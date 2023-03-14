# Load packages
library('mvtnorm')
library('rdrobust')
library('parallel')
library('haven')

triangular <- function(x) {
  return((1-abs(x))*(abs(x)<=1))
}

compare_correlation <- function(Z, Y, threshold) {
  correlation_coefficents <- c()
  for (column in 1:dim(Z)[2]) {
    diff_Z_mean <- Z[,column]-mean(Z[,column])
    diff_Y_mean <- Y-mean(Y)
    correlation_coefficents <- append(correlation_coefficents, sum(diff_Z_mean*diff_Y_mean)/sqrt(sum(diff_Z_mean^2)*sum(diff_Y_mean^2)))
  }
  return(which(abs(correlation_coefficents)>=threshold))
}

calculate_correlation_thresholds <- function(Z, Y) {
  threshold_vector <- c()
  for (column in 1:dim(Z)[2]) {
    sigma <- mean(Z[,column]^2*Y^2)-(mean(Z[,column])^2)*(mean(Y)^2)
    sdZ_sdY <- sqrt(var(Z[,column])*var(Y))
    threshold_vector <- append(threshold_vector, (3*sqrt(sigma))/(sqrt(length(indices))*sdZ_sdY))
  }
  return(threshold_vector)
}

start_time <- Sys.time()

## Load Data
# data <- read_dta("./data/Card_analysis_dataset.dta")

## Restrict to sub-sample of workers with a positive unemployment spell
data <- data[data$dunempl5>0,]
attach(data)

## Assign to those people who had no job before, a duration of the previous job of zero
ind <- which(last_job==0)
last_duration[ind] <- 0

# Basic Model
Xpaper_basis=cbind(female,married,austrian,bluecollar,age, age2,lwage,lwage2,
                   endmo_dum2,endmo_dum3,endmo_dum4, endmo_dum5,endmo_dum6,endmo_dum7,
                   endmo_dum8,endmo_dum9,endmo_dum10, endmo_dum11,endmo_dum12,
                   endy_dum3,endy_dum4,endy_dum5, endy_dum6,endy_dum7,endy_dum8,endy_dum9,
                   endy_dum10,endy_dum11, endy_dum12,endy_dum13,endy_dum14,endy_dum15, 
                   endy_dum16, endy_dum17,endy_dum18,endy_dum19, endy_dum20,endy_dum21)

## Full Model
Xpaper_extended <- cbind(Xpaper_basis,firms,experience,exper2,
                         last_job,last_bluec,last_posnonedur,
                         last_recall,last_noneduration,last_breaks,high_ed,
                         iagrmining,icarsales,ihotel,imanufact,iservice,itransport,iwholesale,
                         reg_dum2, reg_dum3, reg_dum4, reg_dum5,reg_dum6)

## Outcome
Y <- wage_change

## Running Variable
X <- dten1

## Indices of all observations with no entry of NA
indices <- which(apply(!is.na(cbind(X, Y, Xpaper_extended)), 1, all))


# Number of examinations
number_of_examinations <- 8
# Library to use for RDD
# robust - RDRobust
# honest - RDHonest
rdd_library <- "honest"
# Type of estimator (only relevant for RDRobust)
# 1 - Conventional
# 2 - Bias-Corrected
# 3 - Robust
estimator_type <- 1

perform_rdd <- function(X, Y, Z_basis, Z_extended) {
  T = 0+(X >= 0)
  
  coef_vector <- array(NA, c(number_of_examinations))
  standard_error_vector <- array(NA, c(number_of_examinations))
  ci_lower <- array(NA, c(number_of_examinations))
  ci_upper <- array(NA, c(number_of_examinations))
  number_of_covariates <- array(NA, c(number_of_examinations))
  
  n <- length(Y)
  
  covs_with_correlation_greater_002 <- compare_correlation(Z_extended, Y, 0.02)
  covs_with_correlation_greater_005 <- compare_correlation(Z_extended, Y, 0.05)
  covs_with_correlation_greater_01 <- compare_correlation(Z_extended, Y, 0.1)
  covs_with_correlation_greater_02 <- compare_correlation(Z_extended, Y, 0.2)
  covs_with_correlation_greater_calculated_threshold <- compare_correlation(Z_extended, Y, calculate_correlation_thresholds(Z_extended, Y))
  Z_002 <- if (length(covs_with_correlation_greater_002) == 0) NA else Z_extended[,covs_with_correlation_greater_002]
  Z_005 <- if (length(covs_with_correlation_greater_005) == 0) NA else Z_extended[,covs_with_correlation_greater_005]
  Z_01 <- if (length(covs_with_correlation_greater_01) == 0) NA else Z_extended[,covs_with_correlation_greater_01]
  Z_02 <- if (length(covs_with_correlation_greater_02) == 0) NA else Z_extended[,covs_with_correlation_greater_02]
  Z_calculated_threshold <- if (length(covs_with_correlation_greater_calculated_threshold) == 0) NA else Z_extended[,covs_with_correlation_greater_calculated_threshold]
  
  covariate_settings <- list(Z_calculated_threshold, Z_02, Z_01, Z_005, Z_002, NA, Z_basis, Z_extended)
  
  if (rdd_library == "robust") {
    counter <- 1
    for (covariates in covariate_settings) {
      if (all(!is.na(covariates))) {
        rd <- rdrobust(Y, X, covs = covariates)
        number_of_covariates[counter] = dim(covariates)[2]
      } else {
        rd <- rdrobust(Y, X)
        number_of_covariates[counter] = 0
      }
      standard_error_vector[counter] <- rd$se[estimator_type]
      coef_vector[counter] <- rd$coef[estimator_type]
      ci_lower[counter] <- rd$ci[estimator_type,1]
      ci_upper[counter] <- rd$ci[estimator_type,2]
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
        number_of_covariates[counter] = dim(covariates)[2]
      } else {
        Ytilde = Y
        number_of_covariates[counter] = 0
      }
      rd <- RDHonest::RDHonest(Ytilde~X)
      
      standard_error_vector[counter] <- rd$coefficients$std.error
      coef_vector[counter] <- rd$coefficients$estimate
      ci_lower[counter] <- rd$coefficients$conf.low
      ci_upper[counter] <- rd$coefficients$conf.high
      counter = counter + 1
    }
  }
  
  return(matrix(c(number_of_covariates ,coef_vector, standard_error_vector, ci_lower, ci_upper), number_of_examinations, 5))
}

# Results
results <- perform_rdd(X[indices], Y[indices], Xpaper_basis[indices,], Xpaper_extended[indices,])
dimnames(results) <- list(list("Cor.>calc.Th.", "Cor.>0.2", "Cor.>0.1", "Cor.>0.05", "Cor.>0.02", "0 Covs", "Basic Covs", "Extended Covs"), list("#Covs", "Estimator", "Avg. SE", "CI Lower", "CI Upper"))

end_time <- Sys.time()
ellapsed_time <- end_time - start_time
cat("Ellapsed time:", ellapsed_time)

