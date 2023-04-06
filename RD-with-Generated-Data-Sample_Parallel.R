###############################################################################
# This file generates the result of Section 9.2.1 and 9.2.2.                  #
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
#  - The packages mvtnorm, rdrobust and parallel need to be installed.        #
#    Those only need to be installed manually in case the below installation  #
#    does not work properly due to incompatibilities.                         #
###############################################################################


###############################################################################
# Set working directory as well as install and load packages / dependencies   #
###############################################################################
setwd(getSrcDirectory(function(){})[1])

list.of.packages <- c("mvtnorm", "rdrobust", "parallel", "utils")
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

source('R/functions.R')
source('R/RDD_functions.R')


###############################################################################
# Settings for simulation                                                     #
###############################################################################

# Configure number of Monte Carlo replications
# Setting for Section 9.2.1: 10000
# Setting for Section 9.2.2: 10000
number_of_montecarlo_replications <- 100

# Configure sample size of the generated data
# Setting for Section 9.2.1: 1000
# Setting for Section 9.2.2: 10000
n <-1000

# Configure which library to use for RD analysis
#   honest - for RDHonest
#   robust - for RDRobust
# Setting for Section 9.2.1: both were used
# Setting for Section 9.2.2: "honest"
rdd_library <- "honest"

# Configure the type of estimator. This is only relevant if RDRobust is used.
#   1 - Conventional
#   2 - Bias-Corrected
#   3 - Robust
# Setting for Section 9.2.1: 1
# Setting for Section 9.2.2: RDRobust was not used
rdrobust_estimator_type <- 1


###############################################################################
# Executing the simulation                                                    #
###############################################################################

# Parallel execution of the Monte Carlo replications
# A different functions needs to be used dependent on the operating system.
if (Sys.info()['sysname'] == "Windows") {
  # Use number of available kernels except for one
  cluster <- makeCluster(detectCores()-1)
  # Load required packages into the cluster
  load_packages <- parLapply(cluster, 1:length(cluster),
                             function(run) {
                               require('rdrobust')
                               require('mvtnorm')
                             })
  # Export necessary functions to the cluster
  clusterExport(cluster, c('perform_rdd', 'triangular', 'calculate_correlation',
                           'remove_covs_calculated_threshold', 'compare_correlation',
                           'calculate_correlation_thresholds', 'calculate_correlation_threshold_matrix',
                           'remove_covs_with_high_correlation'))
  # Start the parallel simulation
  results_list_of_matrices <- parLapply(cluster, 1:number_of_montecarlo_replications, perform_rdd,
                                        sample_size = n,
                                        rdd_library = rdd_library,
                                        estimator_type = rdrobust_estimator_type)
  stopCluster(cluster)
} else {
  # Start the parallel simulation on all available kernels except for one
  results_list_of_matrices <- mclapply(1:number_of_montecarlo_replications, perform_rdd,
                                       sample_size = n,
                                       rdd_library = rdd_library,
                                       estimator_type = rdrobust_estimator_type,
                                       mc.cores = detectCores()-1)
}

# Transform and split the results
results_list <- lapply(results_list_of_matrices, function(x) { x <- x[[1]] })
selection_list <- lapply(results_list_of_matrices, function(x) { x <- x[[2]] })
results_matrix <- simplify2array(results_list)
selection_matrix <- simplify2array(selection_list)

# Variable which stores the number of different settings of covariate selections
number_of_examinations <- 13

# Initialize vectors for result storage
number_of_covs <- c()
mean_of_estimator <- c()
bias <- c()
standard_deviation <- c()
standard_error <- c()
ci_length <- c()
coverage <- c()
results <- matrix(NA, number_of_examinations, 6, dimnames = list(list("Cor.>calc.Th.", "Th+Del", "Th+Del Simple", "Cor.>0.2", "Cor.>0.1", "Cor.>0.05", "Cor.>0.02", "0 Covs", "1 Cov", "10 Covs", "30 Covs", "50 Covs", "Opt. Cov"), list("#Covs", "Bias", "SD", "Avg. SE", "CI Length", "Coverage")))

# Store results in a matrix
# Iterate over each of the covariate settings examined
for (l in 1:number_of_examinations) {
  # Average the number of selected covariates over all executions
  number_of_covs <- append(number_of_covs, mean(results_matrix[l,5,]))
  # Average the estimation of the average treatment effect over all executions
  mean_of_estimator <- append(mean_of_estimator, mean(results_matrix[l,2,]))
  # Calculate the bias of the estimator
  bias <- append(bias, mean_of_estimator[l]-0.02)
  # Calculate the standard deviation of the estimator
  standard_deviation <- append(standard_deviation, sqrt(mean((results_matrix[l,2,]-mean_of_estimator[l])^2)))
  # Average the standard error over all executions
  standard_error <- append(standard_error, mean(results_matrix[l,1,]))
  # Average the confidence interval length over all executions
  ci_length <- append(ci_length, mean(results_matrix[l,3,]))
  # Calculate the coverage of the confidence intervals
  coverage <- append(coverage, mean(results_matrix[l,4,])*100)
  
  # Store the results in a matrix structure
  results[l,] <- c(number_of_covs[l], bias[l], standard_deviation[l], standard_error[l], ci_length[l], coverage[l])
}

# Calculate percentages, with which each covariate got chosen by the respective selection procedures
selection <- rowSums(selection_matrix, dims = 2)*100/number_of_montecarlo_replications

# Print results on covariate selection
selection

# Print inference results
results
