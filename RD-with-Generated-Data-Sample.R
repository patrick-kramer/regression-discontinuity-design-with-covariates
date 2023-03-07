# Load packages
library('mvtnorm')
library('rdrobust')

start_time <- Sys.time()

number_of_montecarlo_replications <- 10
# Type of RDRobust estimator
# 1 - Conventional
# 2 - Bias-Corrected
# 3 - Robust
type_of_estimator <- 2

coef_vector <- array(NA, c(number_of_montecarlo_replications, 6))
standard_error_vector <- array(NA, c(number_of_montecarlo_replications, 6))
ci_length_vector <- array(NA, c(number_of_montecarlo_replications, 6))
number_for_coverage <- array(0, c(6))

for (n in 1:number_of_montecarlo_replications) {
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
  
  # RDD without covariates
  rd_without_covs <- rdrobust(Y, X)
  standard_error_vector[n, 1] <- rd_without_covs$se[type_of_estimator]
  coef_vector[n, 1] <- rd_without_covs$coef[type_of_estimator]
  ci_length_vector[n, 1] <- rd_without_covs$ci[type_of_estimator,2]-rd_without_covs$ci[type_of_estimator,1]
  if (0.02>=rd_without_covs$ci[type_of_estimator,1] && 0.02<=rd_without_covs$ci[type_of_estimator,2]) {
    number_for_coverage[1] = number_for_coverage[1] + 1
  }
  
  # RDD with 1 covariate
  rd_one_cov <- rdrobust(Y, X, covs = matrix_Z[,1])
  standard_error_vector[n, 2] <- rd_one_cov$se[type_of_estimator]
  coef_vector[n, 2] <- rd_one_cov$coef[type_of_estimator]
  ci_length_vector[n, 2] <- rd_one_cov$ci[type_of_estimator,2]-rd_one_cov$ci[type_of_estimator,1]
  if (0.02>=rd_one_cov$ci[type_of_estimator,1] && 0.02<=rd_one_cov$ci[type_of_estimator,2]) {
    number_for_coverage[2] = number_for_coverage[2] + 1
  }
  
  # RDD with 10 covariates
  rd_ten_covs <- rdrobust(Y, X, covs = matrix_Z[,1:10])
  standard_error_vector[n, 3] <- rd_ten_covs$se[type_of_estimator]
  coef_vector[n, 3] <- rd_ten_covs$coef[type_of_estimator]
  ci_length_vector[n, 3] <- rd_ten_covs$ci[type_of_estimator,2]-rd_ten_covs$ci[type_of_estimator,1]
  if (0.02>=rd_ten_covs$ci[type_of_estimator,1] && 0.02<=rd_ten_covs$ci[type_of_estimator,2]) {
    number_for_coverage[3] = number_for_coverage[3] + 1
  }
  
  # RDD with 30 covariates
  rd_thirty_covs <- rdrobust(Y, X, covs = matrix_Z[,1:30])
  standard_error_vector[n, 4] <- rd_thirty_covs$se[type_of_estimator]
  coef_vector[n, 4] <- rd_thirty_covs$coef[type_of_estimator]
  ci_length_vector[n, 4] <- rd_thirty_covs$ci[type_of_estimator,2]-rd_thirty_covs$ci[type_of_estimator,1]
  if (0.02>=rd_thirty_covs$ci[type_of_estimator,1] && 0.02<=rd_thirty_covs$ci[type_of_estimator,2]) {
    number_for_coverage[4] = number_for_coverage[4] + 1
  }
  
  # RDD with 50 covariates
  rd_fifty_covs <- rdrobust(Y, X, covs = matrix_Z[,1:50])
  standard_error_vector[n, 5] <- rd_fifty_covs$se[type_of_estimator]
  coef_vector[n, 5] <- rd_fifty_covs$coef[type_of_estimator]
  ci_length_vector[n, 5] <- rd_fifty_covs$ci[type_of_estimator,2]-rd_fifty_covs$ci[type_of_estimator,1]
  if (0.02>=rd_fifty_covs$ci[type_of_estimator,1] && 0.02<=rd_fifty_covs$ci[type_of_estimator,2]) {
    number_for_coverage[5] = number_for_coverage[5] + 1
  }
  
  # RDD with optimal covariate
  rd_optimal_cov <- rdrobust(Y, X, covs = Z_times_alpha)
  standard_error_vector[n, 6] <- rd_optimal_cov$se[type_of_estimator]
  coef_vector[n, 6] <- rd_optimal_cov$coef[type_of_estimator]
  ci_length_vector[n, 6] <- rd_optimal_cov$ci[type_of_estimator,2]-rd_optimal_cov$ci[type_of_estimator,1]
  if (0.02>=rd_optimal_cov$ci[type_of_estimator,1] && 0.02<=rd_optimal_cov$ci[type_of_estimator,2]) {
    number_for_coverage[6] = number_for_coverage[6] + 1
  }
  
  cat("Completed run", n, "\n")
}

mean_of_estimator <- c()
bias <- c()
standard_deviation <- c()
standard_error <- c()
ci_length <- c()
coverage <- c()
results <- matrix(NA, 6, 5, dimnames = list(list("0 Covs", "1 Cov", "10 Covs", "30 Covs", "50 Covs", "Opt. Cov"), list("Bias", "SD", "Avg. SE", "CI Length", "Coverage")))

for (l in 1:6) {
  mean_of_estimator <- append(mean_of_estimator, mean(coef_vector[,l]))
  bias <- append(bias, mean_of_estimator[l]-0.02)
  standard_deviation <- append(standard_deviation, sqrt(mean((coef_vector[,l]-mean_of_estimator[l])^2)))
  standard_error <- append(standard_error, mean(standard_error_vector[,l]))
  ci_length <- append(ci_length, mean(ci_length_vector[,l]))
  coverage <- append(coverage, number_for_coverage[l]*100/number_of_montecarlo_replications)
  
  results[l,] <- c(bias[l], standard_deviation[l], standard_error[l], ci_length[l], coverage[l])
}

results

end_time <- Sys.time()
ellapsed_time <- end_time - start_time
cat("Ellapsed time:", ellapsed_time)
