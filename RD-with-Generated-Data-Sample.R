# Load packages
library('mvtnorm')
library('rdrobust')

# Generate finite sample data
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

Y_0 <- Z[,1]+0.36+0.96*X+5.47*X^2+15.28*X^3+15.87*X^4+5.14*X^5+0.22*Z_times_alpha
Y_1 <- Z[,1]+0.38+0.62*X-2.84*X^2+8.42*X^3-10.24*X^4+4.31*X^5+0.28*Z_times_alpha

Y <- (1-T)*Y_0+T*Y_1

# Perform RDD without Covariates
print("--- RDD without covariates ---")
rd <- rdrobust(Y, X)
rd$bias
rd$se
rd$ci
rd$coef
rd$V_cl_l
rd$V_cl_r

print("--- RDD with 10 most important covariates")
rd_ten_covs <- rdrobust(Y, X, covs = matrix_Z[,1:10])
rd_ten_covs$bias
rd_ten_covs$se
rd_ten_covs$ci
rd_ten_covs$coef
rd_ten_covs$V_cl_l
rd_ten_covs$V_cl_r