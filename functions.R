triangular <- function(x) {
  return((1-abs(x))*(abs(x)<=1))
}

calculate_correlation <- function(A, B) {
  diff_A_mean <- A-mean(A)
  diff_B_mean <- B-mean(B)
  correlation_coefficent <- sum(diff_A_mean*diff_B_mean)/sqrt(sum(diff_A_mean^2)*sum(diff_B_mean^2))
  return(correlation_coefficent)
}

compare_correlation <- function(Z, Y, threshold) {
  correlation_coefficents <- c()
  for (column in 1:dim(Z)[2]) {
    correlation_coefficents <- append(correlation_coefficents, calculate_correlation(Z[,column], Y))
  }
  return(which(abs(correlation_coefficents)>=threshold))
}

calculate_correlation_thresholds <- function(Z, Y, factor, data_size) {
  threshold_vector <- c()
  for (column in 1:dim(Z)[2]) {
    sigma <- mean(Z[,column]^2*Y^2)-(mean(Z[,column])^2)*(mean(Y)^2)
    sdZ_sdY <- sqrt(var(Z[,column])*var(Y))
    threshold_vector <- append(threshold_vector, (factor*sqrt(sigma))/(sqrt(data_size)*sdZ_sdY))
  }
  return(threshold_vector)
}

calculate_correlation_threshold_matrix <- function(Z, factor, data_size) {
  threshold_matrix <- matrix(0, nrow=ncol(Z), ncol = ncol(Z))
  for (i in 1:ncol(Z)) {
    for (j in i:ncol(Z)) {
      sigma <- mean(Z[,i]^2*Z[,j]^2)-(mean(Z[,i])^2)*(mean(Z[,j])^2)
      sdZ <- sqrt(var(Z[,i])*var(Z[,j]))
      corr_coeff <- (factor*sqrt(sigma))/(sqrt(data_size)*sdZ)
      threshold_matrix[i,j] <- corr_coeff
      threshold_matrix[j,i] <- corr_coeff
    }
  }
  return(threshold_matrix)
}

remove_covs_with_high_correlation <- function(Z, number) {
  if (isTRUE(ncol(Z)>=1)) {
    covariates_to_delete <- c()
    cor_matrix <- abs(cor(Z))
    covariates_to_delete <- append(covariates_to_delete, which(colSums(is.na(cor_matrix)) == nrow(cor_matrix)-1))
    cor_matrix[,covariates_to_delete] <- 0
    cor_matrix[covariates_to_delete,] <- 0
    number <- number - length(covariates_to_delete)
    if (number > 0) {
      for (i in 1:nrow(cor_matrix)) {
        for (j in 1:i) {
          cor_matrix[i,j] <- 0
        }
      }
      for (i in 1:number) {
        index <- which(cor_matrix == max(cor_matrix), arr.ind = TRUE)
        covariates_to_delete <- append(covariates_to_delete, index[[1,2]])
        cor_matrix[index[1,2],] <- 0
        cor_matrix[,index[1,2]] <- 0
      }
    }
    warning("Removed ", length(covariates_to_delete), " covariates.")
    return(Z[,-covariates_to_delete])
  } else {
    return(Z)
  }
}

remove_covs_with_correlation_larger_threshold <- function(Z, threshold) {
  if (isTRUE(ncol(Z)>1)) {
    cor_matrix <- cor(Z)
    pos <- which(abs(cor_matrix) >= threshold, arr.ind = TRUE)
    pos_ordered <- pos[order(pos[,1], pos[,2]),]
    indices_to_delete_with_duplicates <- c()
    for (i in 1:nrow(pos_ordered)) {
      if (pos_ordered[i,1] < pos_ordered[i,2] && !(pos_ordered[i,1] %in% indices_to_delete_with_duplicates)) {
        indices_to_delete_with_duplicates <- append(indices_to_delete_with_duplicates, pos_ordered[i,2])
      }
    }
    indices_to_delete <- indices_to_delete_with_duplicates[!duplicated(indices_to_delete_with_duplicates)]
    if (length(indices_to_delete)>0) {
      warning("Covariates with high correlation detected. Deleted ", length(indices_to_delete), "duplicates.")
      return(Z[,-indices_to_delete])
    } else {
      return(Z)
    }
  } else {
    return(Z)
  }
}

remove_covs_calculated_threshold <- function(Z, Y, factor, data_size, simple_deletion = TRUE) {
  if (isTRUE(ncol(Z)>1)) {
    cor_matrix <- cor(Z)
    pos <- which(abs(cor_matrix) >= calculate_correlation_threshold_matrix(Z, factor, data_size), arr.ind = TRUE)
    if (simple_deletion) {
      pos_ordered <- pos
    } else {
      pos_ordered <- pos[order(abs(vapply(pos[,1], function(x) { calculate_correlation(Z[,x],Y) }, numeric(1))), decreasing = TRUE),]
    }
    indices_to_delete_with_duplicates <- c()
    for (i in 1:nrow(pos_ordered)) {
      if (simple_deletion) {
        if (pos_ordered[i,1] < pos_ordered[i,2]) {
          correlation_one <- abs(calculate_correlation(Z[,pos_ordered[i,1]], Y))
          correlation_two <- abs(calculate_correlation(Z[,pos_ordered[i,2]], Y))
          if (correlation_one <= correlation_two) {
            indices_to_delete_with_duplicates <- append(indices_to_delete_with_duplicates, pos_ordered[i,1])
          } else {
            indices_to_delete_with_duplicates <- append(indices_to_delete_with_duplicates, pos_ordered[i,2])
          }
        }
      } else {
        if ((pos_ordered[i,1] != pos_ordered[i,2]) && !(pos_ordered[i,1] %in% indices_to_delete_with_duplicates)) {
          indices_to_delete_with_duplicates <- append(indices_to_delete_with_duplicates, pos_ordered[i,2])
        }
      }
    }
    indices_to_delete <- indices_to_delete_with_duplicates[!duplicated(indices_to_delete_with_duplicates)]
    if (length(indices_to_delete)>0) {
      warning("Covariates with high correlation detected. Deleted ", length(indices_to_delete), "duplicates.")
      return(Z[,-indices_to_delete])
    } else {
      return(Z)
    }
  } else {
    return(Z)
  }
}

remove_linearly_dependent_covs <- function(Z) {
  Z <- as.matrix(Z)
  ncovs <- ncol(Z)
  if (ncovs >= 2) {
    df <- data.frame(z = Z)
    constant <- rep(1,nrow(df))
    tmp <- lm(constant ~ ., data=df)
    to_keep <- tmp$coefficients[!is.na(tmp$coefficients)]
    ncovs_keep <- length(to_keep)
    to_keep <- names(to_keep[2:ncovs_keep])
    return(as.matrix(df[to_keep]))
  } else {
    return(Z)
  }
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

cross_interactions <- function(Z1,Z2,ident="") {
  p1 <- dim(Z1)[2]
  p2 <- dim(Z2)[2]
  n <- dim(Z1)[1]
  
  out <- matrix(NA,nrow=n,ncol=p1*p2)
  colnames(out) <- rep("a",p1*p2)
  col_count <- 1
  for(i in 1:p1) {
    for(j in 1:p2) {
      out[,col_count] <- Z1[,i]*Z2[,j]
      colnames(out)[col_count] <- sprintf("CI %s: %d * %d",ident,i,j)
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
