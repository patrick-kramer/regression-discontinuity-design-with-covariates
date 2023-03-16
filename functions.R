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

calculate_correlation_thresholds <- function(Z, Y, data_size) {
  threshold_vector <- c()
  for (column in 1:dim(Z)[2]) {
    sigma <- mean(Z[,column]^2*Y^2)-(mean(Z[,column])^2)*(mean(Y)^2)
    sdZ_sdY <- sqrt(var(Z[,column])*var(Y))
    threshold_vector <- append(threshold_vector, (2*sqrt(sigma))/(sqrt(data_size)*sdZ_sdY))
  }
  return(threshold_vector)
}

remove_covs_with_high_correlation <- function(Z) {
  if (isTRUE(ncol(Z)>1)) {
    cor_matrix <- cor(Z)
    positions <- which(cor_matrix >= 0.9, arr.ind = TRUE)
    indices_to_delete_with_duplicates <- c()
    for (i in 1:nrow(positions)) {
      if (positions[i,1] < positions[i,2]) {
        indices_to_delete_with_duplicates <- append(indices_to_delete_with_duplicates, positions[i,2])
      }
    }
    indices_to_delete <- indices_to_delete_with_duplicates[!duplicated(indices_to_delete_with_duplicates)]
    warning("Covariates with high correlation detected. Deleted ", length(indices_to_delete), "duplicates.")
    return(Z[,-indices_to_delete])
  } else {
    return(Z)
  }
}

remove_linearly_dependent_covs <- function(Z) {
  Z <- as.matrix(Z)
  ncovs <- ncol(Z)
  df <- data.frame(z = Z)
  constant <- rep(1,nrow(df))
  tmp <- lm(constant ~ ., data=df)
  to_keep <- tmp$coefficients[!is.na(tmp$coefficients)]
  ncovs_keep <- length(to_keep)
  to_keep <- names(to_keep[2:ncovs_keep])
  return(list(covs=as.matrix(df[to_keep]), ncovs=ncovs_keep-1))
}

interaction_terms <- function(Z) {
  p <- dim(Z)[2]
  if (p==1) { stop("Input is one-dimensional, no interaction terms computed.") }
  
  out <- matrix(NA,nrow=dim(Z)[1],ncol=p*(p-1)/2)
  colnames(out) <- rep("a",p*(p-1)/2)
  col_count <- 1
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      out[,col_count] <- Z[,i]*Z[,j]
      colnames(out)[col_count] <- sprintf("IT %d*%d",i,j)
      col_count <- col_count+1
    }
  }
  return(out)
}

fourier_basis <- function(Z, order) {
  p <- dim(Z)[2]
  
  out <- matrix(NA,nrow=dim(Z)[1],ncol=p*2*order)
  
  colnames(out) <- rep("a",p*2*order)
  col_count <- 1
  for (i in 1:p) {
    for (j in 1:order) {
      out[,col_count] <- sin(2*pi*j*Z[,i])
      colnames(out)[col_count] <- sprintf("FB %d sin %d",i,j)
      col_count <- col_count+1
      out[,col_count] <- cos(2*pi*j*Z[,i])
      colnames(out)[col_count] <- sprintf("FB %d cos %d",i,j)
      col_count <- col_count+1
    }
  }
  return(out)
}
