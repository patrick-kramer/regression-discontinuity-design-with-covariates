# Load packages
library('mvtnorm')
library('rdrobust')
library('parallel')
setwd(getSrcDirectory(function(){})[1])
source('R/functions.R')
source('R/RDD_functions.R')


#########################################
# Run Simulation                        #
#########################################

start_time <- Sys.time()

# Number of replications
number_of_montecarlo_replications <- 100

# Sample size n
n <-1000

# Library to use for RDD
# robust - RDRobust
# honest - RDHonest
rdd_library <- "robust"

# Type of estimator (only relevant for RDRobust)
# 1 - Conventional
# 2 - Bias-Corrected
# 3 - Robust
rdrobust_estimator_type <- 1

# Number of examinations
number_of_examinations <- 13


# Perform parallel RDD
if (Sys.info()['sysname'] == "Windows") {
  cluster <- makeCluster(11)
  load_packages <- parLapply(cluster, 1:length(cluster),
                             function(run) {
                               require('rdrobust')
                               require('mvtnorm')
                             })
  clusterExport(cluster, c('perform_rdd', 'triangular', 'calculate_correlation',
                           'remove_covs_calculated_threshold', 'compare_correlation',
                           'calculate_correlation_thresholds', 'calculate_correlation_threshold_matrix',
                           'remove_covs_with_high_correlation'))
  results_list_of_matrices <- parLapply(cluster, 1:number_of_montecarlo_replications, perform_rdd,
                                        sample_size = n,
                                        rdd_library = rdd_library,
                                        estimator_type = rdrobust_estimator_type)
  stopCluster(cluster)
} else {
  results_list_of_matrices <- mclapply(1:number_of_montecarlo_replications, perform_rdd,
                                       sample_size = n,
                                       rdd_library = rdd_library,
                                       estimator_type = rdrobust_estimator_type,
                                       mc.cores = 11)
}
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
results <- matrix(NA, number_of_examinations, 6, dimnames = list(list("Cor.>calc.Th.", "Th+Del", "Th+Del Simple", "Cor.>0.2", "Cor.>0.1", "Cor.>0.05", "Cor.>0.02", "0 Covs", "1 Cov", "10 Covs", "30 Covs", "50 Covs", "Opt. Cov"), list("#Covs", "Bias", "SD", "Avg. SE", "CI Length", "Coverage")))

# Store results in a matrix
for (l in 1:number_of_examinations) {
  number_of_covs <- append(number_of_covs, mean(results_matrix[l,5,]))
  mean_of_estimator <- append(mean_of_estimator, mean(results_matrix[l,2,]))
  bias <- append(bias, mean_of_estimator[l]-0.02)
  standard_deviation <- append(standard_deviation, sqrt(mean((results_matrix[l,2,]-mean_of_estimator[l])^2)))
  standard_error <- append(standard_error, mean(results_matrix[l,1,]))
  ci_length <- append(ci_length, mean(results_matrix[l,3,]))
  coverage <- append(coverage, mean(results_matrix[l,4,])*100)
  
  results[l,] <- c(number_of_covs[l], bias[l], standard_deviation[l], standard_error[l], ci_length[l], coverage[l])
}

# Sum up selection matrix
selection <- rowSums(selection_matrix, dims = 2)*100/number_of_montecarlo_replications

# Print results
results

end_time <- Sys.time()
ellapsed_time <- end_time - start_time
cat("Ellapsed time:", ellapsed_time)

