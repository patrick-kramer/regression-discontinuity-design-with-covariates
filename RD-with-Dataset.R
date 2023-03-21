# Load packages
library('mvtnorm')
library('rdrobust')
library('parallel')
library('haven')
setwd(getSrcDirectory(function(){})[1])
source('functions.R')

perform_rdd <- function(X, Y, Z) {
  
  # Define treatment variable
  T = 0+(X >= 0)
  
  # Store number of data entries
  n <- length(Y)
  
  # Initialize vectors for storing results
  coef_vector <- array(NA, c(number_of_examinations))
  standard_error_vector <- array(NA, c(number_of_examinations))
  ci_lower <- array(NA, c(number_of_examinations))
  ci_upper <- array(NA, c(number_of_examinations))
  ci_length <- array(NA, c(number_of_examinations))
  number_of_covariates <- array(NA, c(number_of_examinations))
  
  # Select covariates
  Z_02 <- Z[,compare_correlation(Z, Y, 0.2)]
  Z_calculated_threshold <- Z[,compare_correlation(Z, Y, calculate_correlation_thresholds(Z, Y, length(indices)))]
  
  if (rdd_library == "honest") {
    Z_calculated_threshold <- remove_covs_calculated_threshold(Z_calculated_threshold, length(indices))
    Z_extended <- remove_covs_with_high_correlation(Z[,1:60], 3)
  }
  
  covariate_settings <- list(Z_calculated_threshold, Z_02, NA, Z[,1:38], Z_extended)
  
  gc()
  
  if (rdd_library == "robust") {
    counter <- 1
    for (covariates in covariate_settings) {
      if (isTRUE(ncol(covariates)>0)) {
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
      ci_length[counter] <- ci_upper[counter]-ci_lower[counter]
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
        number_of_covariates[counter] = dim(as.matrix(covariates))[2]
      } else {
        Ytilde = Y
        number_of_covariates[counter] = 0
      }
      rd <- RDHonest::RDHonest(Ytilde~X)
      
      standard_error_vector[counter] <- rd$coefficients$std.error
      coef_vector[counter] <- rd$coefficients$estimate
      ci_lower[counter] <- rd$coefficients$conf.low
      ci_upper[counter] <- rd$coefficients$conf.high
      ci_length[counter] <- ci_upper[counter]-ci_lower[counter]
      counter = counter + 1
    }
  }
  
  return(matrix(c(number_of_covariates ,coef_vector, standard_error_vector, ci_lower, ci_upper, ci_length), number_of_examinations, 6))
}

start_time <- Sys.time()

if (!exists('input_data')) {
  # Load Data
  input_data <- read_dta("./data/Card_analysis_dataset.dta")
  
  # Restrict to sub-sample of workers with a positive unemployment spell
  data <- input_data[input_data$dunempl5>0,]
  attach(data)
}

# Assign to those people who had no job before, a duration of the previous job of zero
ind <- which(last_job==0)
last_duration[ind] <- 0

# Basic Model
Xpaper_basis=cbind(female,married,austrian,bluecollar,age, age2,lwage,lwage2,
                   endmo_dum2,endmo_dum3,endmo_dum4, endmo_dum5,endmo_dum6,endmo_dum7,
                   endmo_dum8,endmo_dum9,endmo_dum10, endmo_dum11,endmo_dum12,
                   endy_dum3,endy_dum4,endy_dum5, endy_dum6,endy_dum7,endy_dum8,endy_dum9,
                   endy_dum10,endy_dum11, endy_dum12,endy_dum13,endy_dum14,endy_dum15, 
                   endy_dum16, endy_dum17,endy_dum18,endy_dum19, endy_dum20,endy_dum21)

# Full Model
Xpaper_extended <- cbind(Xpaper_basis,firms,experience,exper2,
                         last_job,last_bluec,last_posnonedur,
                         last_recall,last_noneduration,last_breaks,high_ed,
                         iagrmining,icarsales,ihotel,imanufact,iservice,itransport,iwholesale,
                         reg_dum2, reg_dum3, reg_dum4, reg_dum5,reg_dum6)

# Outcome
Y <- wage_change

# Running Variable
X <- dten1

# Indices of all observations with no entry of NA
indices <- which(apply(!is.na(cbind(X, Y, Xpaper_extended)), 1, all))
indices <- sample(indices, 10000)

rm("data", "input_data")
gc()

# # Calculate interaction terms for basic and extended covariates
# Z_interaction <- interaction_terms(Xpaper_basis[indices,])
# Z_interaction <- Z_interaction[, colSums(Z_interaction != 0) > 0]
# Z_interaction_df <- as.data.frame(Z_interaction)
# Z_interaction_df <- Z_interaction_df[!duplicated(as.list(Z_interaction_df))]
# Z_interaction <- as.matrix(Z_interaction_df)
# 
# # Calculate fourier bases for basic and extended covariates
# Z_fourier <- fourier_basis(Xpaper_extended[indices, c("lwage", "lwage2")], 5)
# 
# # Join the interaction terms and fourier bases
# Z <- cbind(Xpaper_extended[indices,], Z_interaction, Z_fourier)

Z_interaction <- cbind(interaction_terms(cbind(age,age2,lwage,lwage2,experience,exper2,last_noneduration,last_breaks,
                                         firms,female,married,austrian,bluecollar,last_job,last_bluec,last_posnonedur,last_recall,high_ed))[indices,],
                       cross_interactions(cbind(age,age2,lwage,lwage2,experience,exper2,last_noneduration,last_breaks,
                         firms,female,married,austrian,bluecollar,last_job,last_bluec,last_posnonedur,last_recall,high_ed),cbind(endmo_dum2,endmo_dum3,endmo_dum4, endmo_dum5,endmo_dum6,endmo_dum7,
                                                                                                                                 endmo_dum8,endmo_dum9,endmo_dum10, endmo_dum11,endmo_dum12,
                                                                                                                                 endy_dum3,endy_dum4,endy_dum5, endy_dum6,endy_dum7,endy_dum8,endy_dum9,
                                                                                                                                 endy_dum10,endy_dum11, endy_dum12,endy_dum13,endy_dum14,endy_dum15,
                                                                                                                                 endy_dum16, endy_dum17,endy_dum18,endy_dum19, endy_dum20,endy_dum21,
                                                                                                                                 iagrmining,icarsales,ihotel,imanufact,iservice,itransport,iwholesale,
                                                                                                                                 reg_dum2, reg_dum3, reg_dum4, reg_dum5,reg_dum6,
                                                                                                                                 firms2,firms3,firms4, firms5,firms6),ident="Basic")[indices,],
                       cross_interactions(cbind(endmo_dum2,endmo_dum3,endmo_dum4, endmo_dum5,endmo_dum6,endmo_dum7,
                         endmo_dum8,endmo_dum9,endmo_dum10, endmo_dum11,endmo_dum12),cbind(endy_dum3,endy_dum4,endy_dum5, endy_dum6,endy_dum7,endy_dum8,endy_dum9,
                                                                                           endy_dum10,endy_dum11, endy_dum12,endy_dum13,endy_dum14,endy_dum15,
                                                                                           endy_dum16, endy_dum17,endy_dum18,endy_dum19, endy_dum20,endy_dum21,
                                                                                           iagrmining,icarsales,ihotel,imanufact,iservice,itransport,iwholesale,
                                                                                           reg_dum2, reg_dum3, reg_dum4, reg_dum5,reg_dum6,
                                                                                           firms2,firms3,firms4, firms5,firms6),ident="endmo")[indices,],
                        cross_interactions(cbind(endy_dum3,endy_dum4,endy_dum5, endy_dum6,endy_dum7,endy_dum8,endy_dum9,
                         endy_dum10,endy_dum11, endy_dum12,endy_dum13,endy_dum14,endy_dum15,
                         endy_dum16, endy_dum17,endy_dum18,endy_dum19, endy_dum20,endy_dum21),cbind(iagrmining,icarsales,ihotel,imanufact,iservice,itransport,iwholesale,
                                                                                                    reg_dum2, reg_dum3, reg_dum4, reg_dum5,reg_dum6,
                                                                                                    firms2,firms3,firms4, firms5,firms6),ident="endy")[indices,],
                        cross_interactions(cbind(iagrmining,icarsales,ihotel,imanufact,iservice,itransport,iwholesale),cbind(reg_dum2, reg_dum3, reg_dum4, reg_dum5,reg_dum6,
                                                                                                     firms2,firms3,firms4, firms5,firms6),ident="sector")[indices,],
                        cross_interactions(cbind(reg_dum2, reg_dum3, reg_dum4, reg_dum5,reg_dum6),cbind(firms2,firms3,firms4, firms5,firms6),ident="region")[indices,])

Z_interaction <- Z_interaction[, colSums(Z_interaction != 0) > 0]

Z_fourier <- fourier_basis(cbind(lwage,lwage2),5)[indices,]

# Create High-Dimensional Covariate set
Z <- cbind(Xpaper_extended[indices,],firms2[indices],firms3[indices],firms4[indices], firms5[indices],firms6[indices],
           Z_fourier, Z_interaction)

rm('Xpaper_basis', "Xpaper_extended")
gc()

# Number of examinations
number_of_examinations <- 5

# Library to use for RDD
# robust - RDRobust
# honest - RDHonest
rdd_library <- "honest"

# Type of estimator (only relevant for RDRobust)
# 1 - Conventional
# 2 - Bias-Corrected
# 3 - Robust
estimator_type <- 1

# Results
results <- perform_rdd(X[indices], Y[indices], Z)
dimnames(results) <- list(list("Cor.>calc.Th.", "Cor.>0.2", "0 Covs", "Basic Covs", "Extended Covs"), list("#Covs", "Estimator", "Avg. SE", "CI Lower", "CI Upper", "CI Length"))

# Print results
results

end_time <- Sys.time()
ellapsed_time <- end_time - start_time
cat("Ellapsed time:", ellapsed_time)

