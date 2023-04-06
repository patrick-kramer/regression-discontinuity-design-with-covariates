#' This function executes a RD analysis on generated data including covariates
#' chosen by different selection procedures.
#'
#' @param run Indicates the number of execution when being executed parallel
#' @param sample_size The sample size
#' @param rdd_library The R package to use for the RD analysis. Possible values
#'                    are "honest" for the package RDHonest and "robust" for
#'                    the package RDRobust.
#' @param estimator_type This parameter is just relevant when using RDRobust
#'                       (otherwise it can be ignored). It indicates the estimator
#'                       type used in the RD analysis. Possible values are:
#'                       1 - for conventional estimator
#'                       2 - for bias-corrected estimator
#'                       3 - for robust estimator
#'
#' @return This functions returns the results of the RD analysis (bias, standard
#'         deviation, standard error, confidence intervals, coverage) as well as
#'         the results on the selection of covariates
#' @export
perform_rdd <- function(run, sample_size, rdd_library = "honest", estimator_type = 1) {
  
  # Initialize vectors for storing results
  coef_vector <- array(NA, c(13))
  standard_error_vector <- array(NA, c(13))
  ci_length_vector <- array(NA, c(13))
  coverage_vector <- array(0, c(13))
  number_of_covs_vector <- array(NA, c(13))
  selected_covs_vector <- matrix(0, 200, 3)
  
  # Generate finite sample data (score and covariates)
  X <- 2*rbeta(sample_size, 2, 4)-1
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
  Z_epsilon <- rmvnorm(sample_size, sigma = variance_matrix)
  alpha <- 2/(c(1:200))^2
  matrix_Z <- matrix(Z_epsilon[,2:201], sample_size, byrow = FALSE)
  vector_alpha <- matrix(alpha, 200, 1, byrow = TRUE)
  Z_times_alpha <- matrix_Z %*% vector_alpha
  
  # Define potential outcomes
  Y_0 <- Z_epsilon[,1]+0.36+0.96*X+5.47*X^2+15.28*X^3+15.87*X^4+5.14*X^5+0.22*Z_times_alpha
  Y_1 <- Z_epsilon[,1]+0.38+0.62*X-2.84*X^2+8.42*X^3-10.24*X^4+4.31*X^5+0.28*Z_times_alpha
  
  # Define outcome
  Y <- (1-T)*Y_0+T*Y_1
  
  # Set column names of Z
  colnames(matrix_Z) <- c(1:200)
  
  # Select covariates according to different selection procedures
  # Procedure: Correlation threshold of 0.02
  Z_002 <- matrix_Z[,compare_correlation(matrix_Z, Y, 0.02)]
  # Procedure: Correlation threshold of 0.05
  Z_005 <- matrix_Z[,compare_correlation(matrix_Z, Y, 0.05)]
  # Procedure: Correlation threshold of 0.1
  Z_01 <- matrix_Z[,compare_correlation(matrix_Z, Y, 0.1)]
  # Procedure: Correlation threshold of 0.2
  Z_02 <- matrix_Z[,compare_correlation(matrix_Z, Y, 0.2)]
  # Procedure: Calculated correlation threshold (Section 8.3.1 of the thesis)
  selected_indices_threshold_method <- compare_correlation(matrix_Z, Y, calculate_correlation_thresholds(matrix_Z, Y, sample_size))
  Z_calculated_threshold <- matrix_Z[, selected_indices_threshold_method]
  # Procedure: Calculated correlation threshold with simple deletion (Section 8.3.2 of the thesis)
  Z_calculated_threshold_and_deletion_simple <- remove_covs_calculated_threshold(Z_calculated_threshold, Y, sample_size, simple_deletion = TRUE)
  # Procedure: Calculated correlation threshold with advanced deletion (Section 8.3.3 of the thesis)
  Z_calculated_threshold_and_deletion <- remove_covs_calculated_threshold(Z_calculated_threshold, Y, sample_size, simple_deletion = FALSE)
  
  # Raise variable which counts how often a covariate was selected by the respective procedures
  # Counter for selection via calculated threshold (Section 8.3.1):
  for (index in selected_indices_threshold_method) {
    selected_covs_vector[index, 1] <- selected_covs_vector[index, 1] + 1
  }
  # Counter for selection via calculated threshold and simple deletion (Section 8.3.2)
  indices_after_deletion <- as.numeric(colnames(Z_calculated_threshold_and_deletion))
  for (index in indices_after_deletion) {
    selected_covs_vector[index, 2] <- selected_covs_vector[index, 2] + 1
  }
  # Counter for selection via calculated threshold and simple deletion (Section 8.3.3)
  indices_after_deletion_simple <- as.numeric(colnames(Z_calculated_threshold_and_deletion_simple))
  for (index in indices_after_deletion_simple) {
    selected_covs_vector[index, 3] <- selected_covs_vector[index, 3] + 1
  }
  
  # Remove covariates in some cases to fix invertibility issues
  if (sample_size == 1000) {
    Z_002 <- remove_covs_with_high_correlation(Z_002, 45)
    Z_005 <- remove_covs_with_high_correlation(Z_005, 15)
  }
  
  # Store the different settings of covariate selections in a list
  covariate_settings <- list(as.matrix(Z_calculated_threshold),
                             as.matrix(Z_calculated_threshold_and_deletion),
                             as.matrix(Z_calculated_threshold_and_deletion_simple),
                             as.matrix(Z_02),
                             as.matrix(Z_01),
                             as.matrix(Z_005),
                             as.matrix(Z_002),
                             NA,
                             as.matrix(matrix_Z[,1]),
                             matrix_Z[,1:10],
                             matrix_Z[,1:30],
                             matrix_Z[,1:50],
                             Z_times_alpha)
  
  if (rdd_library == "robust") {
    ### This section performs RDD with the package RDRobust ###
    
    counter <- 1
    # Iterate over all covariate settings
    for (covariates in covariate_settings) {
      if (isTRUE(ncol(covariates)>0)) {
        # If covariates were selected, perform RDD with covariates
        rd <- rdrobust(Y, X, covs = covariates)
      } else {
        # If there are no selected covariates, perform RDD without covariates
        rd <- rdrobust(Y, X)
      }
      # Store standard error of the estimation
      standard_error_vector[counter] <- rd$se[estimator_type]
      # Store estimated average treatment effect
      coef_vector[counter] <- rd$coef[estimator_type]
      # Store length of confidence interval
      ci_length_vector[counter] <- rd$ci[estimator_type,2]-rd$ci[estimator_type,1]
      # If the true average treatment effect is inside the interval, raise the counter
      # Later on, this helps to calculate the coverage
      if (0.02>=rd$ci[estimator_type,1] && 0.02<=rd$ci[estimator_type,2]) {
        coverage_vector[counter] = 1
      }
      # Store the number of selected covariates
      if (length(ncol(covariates)) == 0) {
        number_of_covs_vector[counter] <- 0
      } else {
        number_of_covs_vector[counter] <- ncol(covariates)
      }
      counter <- counter + 1
    }
  } else if (rdd_library == "honest") {
    ### This section performs RDD with the package RDHonest ###
    
    counter <- 1
    # Iterate over all covariate settings
    for (covariates in covariate_settings) {
      if (isTRUE(ncol(covariates)>0)) {
        # Store bandwidth which is taken for RDD without covariates
        # By default RDHonest chooses this bandwidth MSE-optimal
        h <- RDHonest::RDHonest(Y~X)$coefficients$bandwidth
        # Calculate the covariate adjustment according to Definition 7.2
        # in Section 7.1
        kernel_weights <- triangular(X/h)/h
        cov_lm <- lm(covariates~1+X+T+X*T, weights = kernel_weights)
        v <- cov_lm$residuals
        sigma_Z <- t(v)%*%(v*kernel_weights)/sample_size
        sigma_ZY <- colSums(as.matrix(v*as.vector(kernel_weights*Y)))/sample_size
        gamma_n <- solve(sigma_Z)%*%sigma_ZY
        Ytilde <- Y-covariates%*%gamma_n
      } else {
        # If there are no selected covariates, we do not need to adjust the outcome
        Ytilde = Y
      }
      # Perform the RD analysis
      rdd <- RDHonest::RDHonest(Ytilde~X)
      
      # Store standard error
      standard_error_vector[counter] <- rdd$coefficients$std.error
      # Store estimation of average treatment effect
      coef_vector[counter] <- rdd$coefficients$estimate
      # Store length of the confidence interval
      ci_length_vector[counter] <- rdd$coefficients$conf.high-rdd$coefficients$conf.low
      # If the true average treatment effect is inside the interval, raise the counter
      # Later on, this helps to calculate the coverage
      if (0.02>=rdd$coefficients$conf.low && 0.02<=rdd$coefficients$conf.high) {
        coverage_vector[counter] = 1
      }
      # Store the number of selected covariates
      if (length(ncol(covariates)) == 0) {
        number_of_covs_vector[counter] <- 0
      } else {
        number_of_covs_vector[counter] <- ncol(covariates)
      }
      counter = counter + 1
    }
  }
  # Store the inference results in one matrix
  results_of_run <- matrix(c(standard_error_vector,
                             coef_vector,
                             ci_length_vector,
                             coverage_vector,
                             number_of_covs_vector), 13, 5)
  
  # Return the inference results as well as the statistic on selected covariates
  return(list(res = results_of_run, sel = selected_covs_vector))
}