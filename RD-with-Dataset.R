###############################################################################
# This file generates the result of Section 9.4.                              #
# The respective necessary settings are commented below.                      #
###############################################################################
# In order that this file works, the following has to be satisfied:           #
#  - R must be installed                                                      #
#  - The working directory must be set to the folder which contains this file.#
#    This just has to be done manually in case the below command              #
#        setwd(getSrcDirectory(function(){})[1])                              #
#    does not work.                                                           #
#  - The directory where this file is stored must contain a folder 'R'        #
#    which contains the files 'functions.R' and 'RDD_functions.R'             #
#  - The directory where this file is stored must contain a folder 'data'     #
#    which contains the file 'Card_analysis_dataset.dta'                      #
#  - The packages mvtnorm, rdrobust, RDHonest, parallel, haven and utils need #
#    to be installed.                                                         #
#    Those only need to be installed manually in case the below installation  #
#    does not work properly due to incompatibilities.                         #
###############################################################################


###############################################################################
# Set working directory as well as install and load packages / dependencies   #
###############################################################################
list.of.packages <- c("mvtnorm", "rdrobust", "parallel", "utils", "haven")
new.packages <- list.of.packages[!(list.of.packages %in%
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only=TRUE)
if (!("RDHonest" %in% installed.packages())) {
  if (!requireNamespace("remotes")) {
    install.packages("remotes")
  }
  remotes::install_github("kolesarm/RDHonest")
  library("RDHonest")
}

setwd(getSrcDirectory(function(){})[1])
source('R/functions.R')
source('R/RDD_functions.R')


###############################################################################
# Settings for simulation                                                     #
###############################################################################

# Configure which library to use for RD analysis
#   honest - for RDHonest
#   robust - for RDRobust
# Setting for Section 9.4: "honest"
rdd_library <- "honest"

# Configure the type of estimator. This is only relevant if RDRobust is used.
#   1 - Conventional
#   2 - Bias-Corrected
#   3 - Robust
# Setting for Section 9.4: RDRobust was not used
rd_estimator_type <- 1

# Limits the loaded data set to a certain amount of randomly chosen entries.
# This was necessary due to limited RAM capacities.
# Setting for Section 9.4: 100000
data_size_limit <- 10000


###############################################################################
# Load and prepare the data set                                               #
###############################################################################

if (!exists('input_data')) {
  # Load Data
  input_data <- read_dta("data/Card_analysis_dataset.dta")
  
  # Restrict to sample of workers with a positive unemployment spell
  data <- input_data[input_data$dunempl5>0,]
  attach(data)
}

# In case a person did not have a job before, the duration of the previous job is set to zero
ind <- which(last_job==0)
last_duration[ind] <- 0

# Load basic covariates
Xpaper_basis=cbind(female,married,austrian,bluecollar,age, age2,lwage,lwage2,
                   endmo_dum2,endmo_dum3,endmo_dum4, endmo_dum5,endmo_dum6,endmo_dum7,
                   endmo_dum8,endmo_dum9,endmo_dum10, endmo_dum11,endmo_dum12,
                   endy_dum3,endy_dum4,endy_dum5, endy_dum6,endy_dum7,endy_dum8,endy_dum9,
                   endy_dum10,endy_dum11, endy_dum12,endy_dum13,endy_dum14,endy_dum15, 
                   endy_dum16, endy_dum17,endy_dum18,endy_dum19, endy_dum20,endy_dum21)

# Load additional covariates (for extended set)
Xpaper_extended <- cbind(Xpaper_basis,firms,experience,exper2,
                         last_job,last_bluec,last_posnonedur,
                         last_recall,last_noneduration,last_breaks,high_ed,
                         iagrmining,icarsales,ihotel,imanufact,iservice,itransport,iwholesale,
                         reg_dum2, reg_dum3, reg_dum4, reg_dum5,reg_dum6)

# Define outcome
Y <- wage_change

# Define running variable
X <- dten1

# Indices of all observations with no entry of NA
indices <- which(apply(!is.na(cbind(X, Y, Xpaper_extended)), 1, all))
# Limit the data set to a certain number of data entries
indices <- sample(indices, data_size_limit)

# Remove variables that are not needed anymore to free RAM
rm("data", "input_data")
gc()

# Compute non-trivial interaction and cross-interaction terms of the covariates
Z_interaction <- cbind(interaction_terms(cbind(age,age2,lwage,lwage2,experience,exper2,last_noneduration,last_breaks,
                                         firms,female,married,austrian,bluecollar,last_job,last_bluec,last_posnonedur,last_recall,high_ed))[indices,],
                       cross_interactions(cbind(age,age2,lwage,lwage2,experience,exper2,last_noneduration,last_breaks,
                         firms,female,married,austrian,bluecollar,last_job,last_bluec,last_posnonedur,last_recall,high_ed),cbind(endmo_dum2,endmo_dum3,endmo_dum4, endmo_dum5,endmo_dum6,endmo_dum7,
                                                                                                                                 endmo_dum8,endmo_dum9,endmo_dum10, endmo_dum11,endmo_dum12,
                                                                                                                                 endy_dum3,endy_dum4,endy_dum5, endy_dum6,endy_dum7,endy_dum8,endy_dum9,
                                                                                                                                 endy_dum10,endy_dum11, endy_dum12,endy_dum13,endy_dum14,endy_dum15,
                                                                                                                                 endy_dum16, endy_dum17,endy_dum18,endy_dum19, endy_dum20,endy_dum21,
                                                                                                                                 iagrmining,icarsales,ihotel,imanufact,iservice,itransport,iwholesale,
                                                                                                                                 reg_dum2, reg_dum3, reg_dum4, reg_dum5,reg_dum6),ident="Basic")[indices,],
                       cross_interactions(cbind(endmo_dum2,endmo_dum3,endmo_dum4, endmo_dum5,endmo_dum6,endmo_dum7,
                         endmo_dum8,endmo_dum9,endmo_dum10, endmo_dum11,endmo_dum12),cbind(endy_dum3,endy_dum4,endy_dum5, endy_dum6,endy_dum7,endy_dum8,endy_dum9,
                                                                                           endy_dum10,endy_dum11, endy_dum12,endy_dum13,endy_dum14,endy_dum15,
                                                                                           endy_dum16, endy_dum17,endy_dum18,endy_dum19, endy_dum20,endy_dum21,
                                                                                           iagrmining,icarsales,ihotel,imanufact,iservice,itransport,iwholesale,
                                                                                           reg_dum2, reg_dum3, reg_dum4, reg_dum5,reg_dum6),ident="endmo")[indices,],
                        cross_interactions(cbind(endy_dum3,endy_dum4,endy_dum5, endy_dum6,endy_dum7,endy_dum8,endy_dum9,
                         endy_dum10,endy_dum11, endy_dum12,endy_dum13,endy_dum14,endy_dum15,
                         endy_dum16, endy_dum17,endy_dum18,endy_dum19, endy_dum20,endy_dum21),cbind(iagrmining,icarsales,ihotel,imanufact,iservice,itransport,iwholesale,
                                                                                                    reg_dum2, reg_dum3, reg_dum4, reg_dum5,reg_dum6),ident="endy")[indices,],
                        cross_interactions(cbind(iagrmining,icarsales,ihotel,imanufact,iservice,itransport,iwholesale),cbind(reg_dum2, reg_dum3, reg_dum4, reg_dum5,reg_dum6),ident="sector")[indices,])

# Select all covariates that are not equal to the zero vector
Z_interaction <- Z_interaction[, colSums(Z_interaction != 0) > 0]

# Compute trignometric transformations on covariates up to order 5
Z_fourier <- fourier_basis(cbind(lwage,lwage2),5)[indices,]

# Combine all the covariates to one matrix
Z <- cbind(Xpaper_extended[indices,], Z_fourier, Z_interaction)

# Remove variables that are not needed anymore to free RAM
rm('Xpaper_basis', "Xpaper_extended")
gc()


###############################################################################
# Perform RD analysis as well as arrange and print results                    #
###############################################################################

output <- perform_rdd_on_data(X[indices], Y[indices], Z, rdd_library = rdd_library, estimator_type = rd_estimator_type)

# Store results
results <- output[['res']]
dimnames(results) <- list(list("(CCTAD)", "(CCTSD)", "Cor.>0.2", "0 Covs", "Basic Covs", "Extended Covs"), list("#Covs", "Estimator", "Avg. SE", "CI Lower", "CI Upper", "CI Length"))
# Add additional row for the selection procedure (CCT) and store its number of selected covariates
results_with_cct <- rbind(
  matrix(c(output[['number_sel_cov_threshold']], NA, NA, NA, NA, NA), nrow=1, ncol=6, dimnames = list(c("(CCT)"), c())),
  results
)

# Store covariate selection results for the procedures (CCTSD) and (CCTAD)
selection <- output[['sel']]

# Print selection results
selection

# Print inference results
results_with_cct
