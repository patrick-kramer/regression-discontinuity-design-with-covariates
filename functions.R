triangular <- function(x) {
  return((1-abs(x))*(abs(x)<=1))
}

compare_correlation <- function(Z, Y, threshold) {
  correlation_coefficents <- c()
  for (column in 1:200) {
    diff_Z_mean <- Z[,column]-mean(Z[,column])
    diff_Y_mean <- Y-mean(Y)
    correlation_coefficents <- append(correlation_coefficents, sum(diff_Z_mean*diff_Y_mean)/sqrt(sum(diff_Z_mean^2)*sum(diff_Y_mean^2)))
  }
  return(which(abs(correlation_coefficents)>=threshold))
}

calculate_correlation_thresholds <- function(Z, Y, data_size) {
  threshold_vector <- c()
  for (column in 1:200) {
    sigma <- mean(Z[,column]^2*Y^2)-(mean(Z[,column])^2)*(mean(Y)^2)
    sdZ_sdY <- sqrt(var(Z[,column])*var(Y))
    threshold_vector <- append(threshold_vector, (2*sqrt(sigma))/(sqrt(data_size)*sdZ_sdY))
  }
  return(threshold_vector)
}