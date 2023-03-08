# Load packages
library('mvtnorm')
library('rdrobust')
library('parallel')

triangular <- function(x) {
  return((1-abs(x))*(abs(x)<=1))
}

start_time <- Sys.time()

# Number of replications
number_of_montecarlo_replications <- 1000

perform_rdd <- function(n) {
  # Library to use for RDD
  # robust - RDRobust
  # honest - RDHonest
  rdd_library <- "robust"
  
  # Type of estimator (only relevant for RDRobust)
  # 1 - Conventional
  # 2 - Bias-Corrected
  # 3 - Robust
  estimator_type <- 2
  
  coef_vector <- array(NA, c(6))
  standard_error_vector <- array(NA, c(6))
  ci_length_vector <- array(NA, c(6))
  coverage_vector <- array(0, c(6))
  
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
  
  n <- length(Y)
  covariate_settings <- list(NA, matrix_Z[,1], matrix_Z[,1:10], matrix_Z[,1:30], matrix_Z[,1:50], Z_times_alpha)
  
  if (rdd_library == "robust") {
    counter <- 1
    for (covariates in covariate_settings) {
      if (all(!is.na(covariates))) {
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
      counter <- counter + 1
    }
  } else if (rdd_library == "honest") {
    counter <- 1
    for (covariates in covariate_settings) {
      if (all(!is.na(covariates))) {
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
      counter = counter + 1
    }
  }
  
  cat("Completed run", n, "\n")
  
  return(matrix(c(standard_error_vector, coef_vector, ci_length_vector, coverage_vector), 6, 4))
}

results_list_of_matrices <- mclapply(1:number_of_montecarlo_replications, perform_rdd, mc.cores = 10)
results_matrix <- simplify2array(results_list_of_matrices)

mean_of_estimator <- c()
bias <- c()
standard_deviation <- c()
standard_error <- c()
ci_length <- c()
coverage <- c()
results <- matrix(NA, 6, 5, dimnames = list(list("0 Covs", "1 Cov", "10 Covs", "30 Covs", "50 Covs", "Opt. Cov"), list("Bias", "SD", "Avg. SE", "CI Length", "Coverage")))

for (l in 1:6) {
  mean_of_estimator <- append(mean_of_estimator, mean(results_matrix[l,2,]))
  bias <- append(bias, mean_of_estimator[l]-0.02)
  standard_deviation <- append(standard_deviation, sqrt(mean((results_matrix[l,2,]-mean_of_estimator[l])^2)))
  standard_error <- append(standard_error, mean(results_matrix[l,1,]))
  ci_length <- append(ci_length, mean(results_matrix[l,3,]))
  coverage <- append(coverage, mean(results_matrix[l,4,])*100)
  
  results[l,] <- c(bias[l], standard_deviation[l], standard_error[l], ci_length[l], coverage[l])
}

results

end_time <- Sys.time()
ellapsed_time <- end_time - start_time
cat("Ellapsed time:", ellapsed_time)

