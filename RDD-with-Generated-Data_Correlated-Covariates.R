###############################################################################
# This file generates the result of Section 9.3.                              #
# The necessary settings are commented below.                                 #
###############################################################################
# In order that this file works, the following has to be satisfied:           #
#  - R must be installed                                                      #
#  - The working directory must be set to the folder which contains this file.#
#    This just has to be done manually in case the below command              #
#        setwd(getSrcDirectory(function(){})[1])                              #
#    does not work.                                                           #
#  - The directory where this file is stored must contain a folder 'R'        #
#    which contains the files 'functions.R' and 'RDD_functions.R'             #
#  - The packages mvtnorm, rdrobust, RDHonest, parallel and utils need to be  #
#    installed.                                                               #
#    They only need to be installed manually in case the below installation   #
#    does not work properly due to incompatibilities.                         #
###############################################################################


###############################################################################
# Set working directory as well as install and load packages / dependencies   #
###############################################################################
list.of.packages <- c("mvtnorm", "rdrobust", "parallel", "utils", "glmnet")
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

# Configure number of Monte Carlo replications
# Setting for Section 9.3: 10000
number_of_montecarlo_replications <- 100

# Configure sample size of the generated data
# Setting for Section 9.3: 1000
n <-1000

# Configure which library to use for RD analysis
#   honest - for RDHonest
#   robust - for RDRobust
# Setting for Section 9.3: "honest"
rdd_library <- "honest"

# Configure the type of estimator. This is only relevant if RDRobust is used.
#   1 - Conventional
#   2 - Bias-Corrected
#   3 - Robust
# Setting for Section 9.3: RDRobust was not used
rdrobust_estimator_type <- 1


###############################################################################
# Executing the simulation                                                    #
###############################################################################

# Parallel execution of the Monte Carlo replications
# A different functions needs to be used dependent on the operating system
if (Sys.info()['sysname'] == "Windows") {
  # Use number of available kernels except for one
  cluster <- makeCluster(detectCores()-1)
  # Load required packages into the cluster
  load_packages <- parLapply(cluster, 1:length(cluster),
                             function(run) {
                               require('rdrobust')
                               require('mvtnorm')
                               require('glmnet')
                             })
  # Export necessary functions to the cluster
  clusterExport(cluster, c('perform_rdd_redundant_covariates', 'triangular', 'calculate_correlation',
                           'remove_covs_calculated_threshold', 'compare_correlation',
                           'calculate_correlation_thresholds', 'calculate_correlation_threshold_matrix',
                           'remove_covs_with_high_correlation', 'select_covs_via_lasso', 'bstpc', 'BCHtpc'))
  # Start the parallel simulation
  results_list_of_matrices <- parLapply(cluster, 1:number_of_montecarlo_replications, perform_rdd_correlated_covariates,
                                        sample_size = n,
                                        rdd_library = rdd_library,
                                        estimator_type = rdrobust_estimator_type)
  stopCluster(cluster)
} else {
  # Start the parallel simulation on all available kernels except for one
  results_list_of_matrices <- mclapply(1:number_of_montecarlo_replications, perform_rdd_correlated_covariates,
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

# Initialize vectors for result storage
number_of_covs <- c()
mean_of_estimator <- c()
bias <- c()
standard_deviation <- c()
standard_error <- c()
ci_length <- c()
coverage <- c()
results <- matrix(NA, 6, 6, dimnames = list(list("(CCT)", "(CCTAD)", "(CCTSD)", "(Lasso BCH)", "(Lasso OPC)", "(Lasso CV)"),
                                            list("#Covs", "Bias", "SD", "Avg. SE",
                                                 "CI Length", "Coverage")))

# Store results in a matrix
# Iterate over the results of each covariate setting
for (l in 1:6) {
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
  results[l,] <- c(number_of_covs[l], bias[l], standard_deviation[l],
                   standard_error[l], ci_length[l], coverage[l])
}

# Calculate percentages with which each covariate got chosen by the respective
# selection procedures (CCT), (CCTAD) and (CCTSD)
selection <- rowSums(selection_matrix, dims = 2)*100/number_of_montecarlo_replications
# Assign column names to selection matrix
colnames(selection) <- c("(CCT)", "(CCTAD)", "(CCTSD)", "(Lasso BCH)", "(Lasso OPC)", "(Lasso CV)")

# Print results on covariate selection for the procedures (CCT), (CCTAD) and (CCTSD)
selection

# Print inference results
results

