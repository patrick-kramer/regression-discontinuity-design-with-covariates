# Load packages
library('mvtnorm')
library('rdrobust')
library('parallel')
setwd(getSrcDirectory(function(){})[1])
source('functions.R')

#########################################
# Define functions                      #
#########################################

perform_rdd <- function(run) {
  
  # Initialize vectors for storing results
  coef_vector <- array(NA, c(number_of_examinations))
  standard_error_vector <- array(NA, c(number_of_examinations))
  ci_length_vector <- array(NA, c(number_of_examinations))
  coverage_vector <- array(0, c(number_of_examinations))
  number_of_covs_vector <- array(NA, c(number_of_examinations))
  selected_covs_vector <- array(0, c(200))
  
  # Generate finite sample data (score and covariates)
  X <- 2*rbeta(1000, 2, 4)-1
  T <- 0+(X >= 0)
  sigma_z <- 0.1353
  sigma_eps <- 0.1295
  v <- 0.8*sqrt(6)*(sigma_eps)^2/(pi*c(1:200))
  sigma_z_id_matrix <- diag(x = sigma_z^2, 200)
  matrix_rows <- matrix(NA, 200, 201)
  for (k in 1:200) {
    matrix_rows[k,] <- c(v[k], sigma_z_id_matrix[k,])
  }
  variance_matrix <- rbind(c(sigma_eps^2, v), matrix_rows)
  Z <- rmvnorm(1000, sigma = variance_matrix)
  alpha <- 2/(c(1:200))^2
  matrix_Z <- matrix(Z[,2:201], 1000, byrow = FALSE)
  vector_alpha <- matrix(alpha, 200, 1, byrow = TRUE)
  Z_times_alpha <- matrix_Z %*% vector_alpha
  
  # Define potential outcomes
  Y_0 <- Z[,1]+0.36+0.96*X+5.47*X^2+15.28*X^3+15.87*X^4+5.14*X^5+0.22*Z_times_alpha
  Y_1 <- Z[,1]+0.38+0.62*X-2.84*X^2+8.42*X^3-10.24*X^4+4.31*X^5+0.28*Z_times_alpha
  
  # Define outcome
  Y <- (1-T)*Y_0+T*Y_1
  
  # Store number of data entries
  n <- length(Y)

  
  # Select covariates
  Z_002 <- Z[,compare_correlation(matrix_Z, Y, 0.02)]
  Z_005 <- Z[,compare_correlation(matrix_Z, Y, 0.05)]
  Z_01 <- Z[,compare_correlation(matrix_Z, Y, 0.1)]
  Z_02 <- Z[,compare_correlation(matrix_Z, Y, 0.2)]
  selected_indices_threshold_method <- compare_correlation(matrix_Z, Y, calculate_correlation_thresholds(matrix_Z, Y, n))
  Z_calculated_threshold <- Z[, selected_indices_threshold_method]
  
  # Store information about selected covariates
  for (index in selected_indices_threshold_method) {
    selected_covs_vector[index] <- selected_covs_vector[index] + 1
  }
  
  # Remove covariables for invertibility
  Z_002 <- remove_covs_with_high_correlation(Z_002, 40)
  Z_005 <- remove_covs_with_high_correlation(Z_005, 10)
  
  covariate_settings <- list(Z_calculated_threshold, Z_02, Z_01, Z_005, Z_002, NA, as.matrix(matrix_Z[,1]), matrix_Z[,1:10], matrix_Z[,1:30], matrix_Z[,1:50], Z_times_alpha)
  
  if (rdd_library == "robust") {
    counter <- 1
    for (covariates in covariate_settings) {
      if (isTRUE(ncol(covariates)>0)) {
        rd <- rdrobust(Y, X, covs = covariates)
      } else {
        rd <- rdrobust(Y, X)
      }
      standard_error_vector[counter] <- rd$se[estimator_type]
      coef_vector[counter] <- rd$coef[estimator_type]
      ci_length_vector[counter] <- rd$ci[estimator_type,2]-rd$ci[estimator_type,1]
      if (0.02>=rd$ci[estimator_type,1] && 0.02<=rd$ci[estimator_type,2]) {
        coverage_vector[counter] = 1
      }
      if (length(ncol(covariates)) == 0) {
        number_of_covs_vector[counter] <- 0
      } else {
        number_of_covs_vector[counter] <- ncol(covariates)
      }
      counter <- counter + 1
    }
  } else if (rdd_library == "honest") {
    counter <- 1
    for (covariates in covariate_settings) {
      if (isTRUE(ncol(covariates)>0)) {
        h <- RDHonest::RDHonest(Y~X)$coefficients$bandwidth
        kernel_weights <- triangular(X/h)/h
        cov_lm <- lm(covariates~1+X+T+X*T, weights = kernel_weights)
        v <- cov_lm$residuals
        sigma_Z <- t(v)%*%(v*kernel_weights)/n
        sigma_ZY <- colSums(as.matrix(v*as.vector(kernel_weights*Y)))/n
        gamma_n <- solve(sigma_Z)%*%sigma_ZY
        Ytilde <- Y-covariates%*%gamma_n
      } else {
        Ytilde = Y
      }
      rd_ten_covs <- RDHonest::RDHonest(Ytilde~X)
      
      standard_error_vector[counter] <- rd_ten_covs$coefficients$std.error
      coef_vector[counter] <- rd_ten_covs$coefficients$estimate
      ci_length_vector[counter] <- rd_ten_covs$coefficients$conf.high-rd_ten_covs$coefficients$conf.low
      if (0.02>=rd_ten_covs$coefficients$conf.low && 0.02<=rd_ten_covs$coefficients$conf.high) {
        coverage_vector[counter] = 1
      }
      if (length(ncol(covariates)) == 0) {
        number_of_covs_vector[counter] <- 0
      } else {
        number_of_covs_vector[counter] <- ncol(covariates)
      }
      counter = counter + 1
    }
  }
  results_of_run <- matrix(c(standard_error_vector, coef_vector, ci_length_vector, coverage_vector, number_of_covs_vector), number_of_examinations, 5)
  
  return(list(res = results_of_run, sel = selected_covs_vector))
}

#########################################
# Run Simulation                        #
#########################################

start_time <- Sys.time()

# Number of replications
number_of_montecarlo_replications <- 10000

# Library to use for RDD
# robust - RDRobust
# honest - RDHonest
rdd_library <- "honest"

# Type of estimator (only relevant for RDRobust)
# 1 - Conventional
# 2 - Bias-Corrected
# 3 - Robust
estimator_type <- 1

# Number of examinations
number_of_examinations <- 11


# Perform parallized RDD
if (Sys.info()['sysname'] == "Windows") {
  cluster <- makeCluster(11)
  load_packages <- parLapply(cluster, 1:length(cluster),
                             function(run) {
                               require('rdrobust')
                               require('mvtnorm')
                             })
  clusterExport(cluster, c('perform_rdd', 'triangular', 'compare_correlation', 'calculate_correlation_thresholds', 'number_of_examinations'))
  results_list_of_matrices <- parLapply(cluster, 1:number_of_montecarlo_replications, perform_rdd)
  stopCluster(cluster)
} else {
  results_list_of_matrices <- mclapply(1:number_of_montecarlo_replications, perform_rdd, mc.cores = 11)
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
results <- matrix(NA, number_of_examinations, 6, dimnames = list(list("Cor.>calc.Th.", "Cor.>0.2", "Cor.>0.1", "Cor.>0.05", "Cor.>0.02", "0 Covs", "1 Cov", "10 Covs", "30 Covs", "50 Covs", "Opt. Cov"), list("#Covs", "Bias", "SD", "Avg. SE", "CI Length", "Coverage")))

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
selection <- rowSums(selection_matrix)*100/number_of_montecarlo_replications

# Print results
results

end_time <- Sys.time()
ellapsed_time <- end_time - start_time
cat("Ellapsed time:", ellapsed_time)

