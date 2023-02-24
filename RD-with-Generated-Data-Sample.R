# Load packages
library('mvtnorm')
library('rdrobust')

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

# Perform RDD
print("--- RDD without covariates ---")
rd_without_covs <- rdrobust(Y, X)

print("--- RDD with most important covariates ---")
rd_one_cov <- rdrobust(Y, X, covs = matrix_Z[,1])

print("--- RDD with 10 most important covariates ---")
rd_ten_covs <- rdrobust(Y, X, covs = matrix_Z[,1:10])

print("--- RDD with 20 most important covariates ---")
rd_twenty_covs <- rdrobust(Y, X, covs = matrix_Z[,1:20])

print("--- RDD with 50 most important covariates ---")
rd_fifty_covs <- rdrobust(Y, X, covs = matrix_Z[,1:50])

# Collect and print results
# TODO: Add SD and coverage
results <- matrix(NA, 5, 3, dimnames = list(list("0 Covs", "1 Cov", "10 Covs", "20 Covs", "50 Covs"), list("Bias", "Avg. SE", "CI Length")))
results[1,] <- c((rd_without_covs$bias[1]+rd_without_covs$bias[2])/2, rd_without_covs$se[1], rd_without_covs$ci[1,2]-rd_without_covs$ci[1,1])
results[2,] <- c((rd_one_cov$bias[1]+rd_one_cov$bias[2])/2, rd_one_cov$se[1], rd_one_cov$ci[1,2]-rd_one_cov$ci[1,1])
results[3,] <- c((rd_ten_covs$bias[1]+rd_ten_covs$bias[2])/2, rd_ten_covs$se[1], rd_ten_covs$ci[1,2]-rd_ten_covs$ci[1,1])
results[4,] <- c((rd_twenty_covs$bias[1]+rd_twenty_covs$bias[2])/2, rd_twenty_covs$se[1], rd_twenty_covs$ci[1,2]-rd_twenty_covs$ci[1,1])
results[5,] <- c((rd_fifty_covs$bias[1]+rd_fifty_covs$bias[2])/2, rd_fifty_covs$se[1], rd_fifty_covs$ci[1,2]-rd_fifty_covs$ci[1,1])
results
