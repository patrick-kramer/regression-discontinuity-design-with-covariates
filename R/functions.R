#' Computes the triangular kernel
#'
#' @param x The argument passed to the triangular kernel function.
#'
#' @return The value of the triangular kernel function.
#' @export
triangular <- function(x) {
  return((1-abs(x))*(abs(x)<=1))
}

#' Calculate the sample correlation coefficient of two vectors
#'
#' @param A Vector of dimension n
#' @param B Vector of dimension n
#'
#' @return Sample correlation coefficient of A and B
#' @export
#'
#' @examples
#' calculate_correlation(c(1,2,3),c(3,2,1))
calculate_correlation <- function(A, B) {
  diff_A_mean <- A-mean(A)
  diff_B_mean <- B-mean(B)
  correlation_coefficent <- sum(diff_A_mean*diff_B_mean)/sqrt(sum(diff_A_mean^2)*sum(diff_B_mean^2))
  return(correlation_coefficent)
}

#' Calculate correlation thresholds for each covariate with the outcome
#' according to Section 8.3.1 of the thesis
#'
#' @param Z The covariate matrix of dimension nxp
#' @param Y The outcome vector of dimension n
#' @param data_size The sample size n
#'
#' @return The threshold vector containing a threshold for each of the
#'         covariates
#' @export
calculate_correlation_thresholds <- function(Z, Y, data_size) {
  # Store the number of covariates
  p <- ncol(Z)
  # Perform Bonferroni correction
  alpha_ind <- 0.05 / p
  factor <- qnorm(1-(alpha_ind / 2))
  
  threshold_vector <- c()
  # Iterate over all covariates
  for (column in 1:dim(Z)[2]) {
    # Calculate threshold according to the formula in Section 8.3.1
    sigma <- mean(Z[,column]^2*Y^2)-(mean(Z[,column])^2)*(mean(Y)^2)
    sdZ_sdY <- sqrt(var(Z[,column])*var(Y))
    threshold_vector <- append(threshold_vector,
                               (factor*sqrt(sigma))/(sqrt(data_size)*sdZ_sdY))
  }
  return(threshold_vector)
}

#' Calculate sample correlation coefficient of each covariate with the outcome
#' and check whether it is greater or equal a certain threshold
#'
#' @param Z The covariates matrix of dimension nxp
#' @param Y The outcome vector of dimension n
#' @param threshold The threshold vector of dimension n
#'
#' @return Returns indices of covariates which have a correlation to the outcome
#'         greater or equal to the according threshold
#' @export
compare_correlation <- function(Z, Y, threshold) {
  correlation_coefficents <- c()
  
  # Iterate over all columns of the covariate matrix Z
  for (column in 1:dim(Z)[2]) {
    # Append the calculated sample correlation coefficient
    correlation_coefficents <- append(correlation_coefficents,
                                      calculate_correlation(Z[,column], Y))
  }
  return(which(abs(correlation_coefficents)>=threshold))
}

#' Calculate deletion correlation thresholds for each pair of covariates
#' according to Section 8.3.2 of the thesis
#'
#' @param Z A matrix containing the selected covariates by the selection
#'          procedure of Section 8.3.1 of dimension nxs
#' @param data_size The sample size n
#'
#' @return The matrix containing the thresholds for each pair of covariates
#' @export
calculate_correlation_threshold_matrix <- function(Z, data_size) {
  # Save number of covariates
  s <- ncol(Z)
  # Calculate number of pairs of covariates (without considering order)
  N <- s*(s-1)*0.5
  # Perform Bonferroni correction
  alpha_ind <- 0.05 / N
  factor <- qnorm(1-(alpha_ind / 2))
  
  threshold_matrix <- matrix(0, nrow=ncol(Z), ncol = ncol(Z))
  
  # Iterate over all pairs of covariates (without considering order)
  for (i in 1:ncol(Z)) {
    for (j in i:ncol(Z)) {
      # Calculate the deletion tresholds according to Section 8.3.2
      sigma <- mean(Z[,i]^2*Z[,j]^2)-(mean(Z[,i])^2)*(mean(Z[,j])^2)
      sdZ <- sqrt(var(Z[,i])*var(Z[,j]))
      corr_coeff <- (factor*sqrt(sigma))/(sqrt(data_size)*sdZ)
      threshold_matrix[i,j] <- corr_coeff
      threshold_matrix[j,i] <- corr_coeff
    }
  }
  return(threshold_matrix)
}

#' Remove selected covariates according to the procedure described in Section
#' 8.3.2 and 8.3.3 (simple and advanced deletion, respectively)
#'
#' @param Z A matrix containing the selected covariates by the Selection
#'          procedure of Section 8.3.1 of dimension nxs
#' @param Y The outcome vector of dimension of n
#' @param data_size The sample size n
#' @param simple_deletion When set to TRUE, simple deletion is applied (Section
#'                        8.3.2), otherwise advanced deletion (Section 8.3.3)
#'
#' @return The remaining set of covariates after the deletion
#' @export
remove_covs_calculated_threshold <- function(Z, Y, data_size, simple_deletion = TRUE) {
  if (isTRUE(ncol(Z)>1)) {
    # Calculate correlation matrix for the covariates
    cor_matrix <- cor(Z)
    
    # Store positions of the pairs having the absolute value of the sample correlation greater than the
    # calculated threshold
    pos <- which(abs(cor_matrix) >= calculate_correlation_threshold_matrix(Z, data_size), arr.ind = TRUE)
    
    if (simple_deletion) {
      # In case of simple deletion, just keep all those positions
      pos_ordered <- pos
    } else {
      # In case of advanced deletion, sort them according to the correlation to the outcome
      pos_ordered <- pos[order(abs(vapply(pos[,1],
                                          function(x) { calculate_correlation(Z[,x],Y) },
                                          numeric(1))),
                               decreasing = TRUE),]
    }
    
    indices_to_delete_with_duplicates <- c()
    
    # Iterate over all stored positions
    for (i in 1:nrow(pos_ordered)) {
      if (simple_deletion) {
        # In case of simple deletion, delete the covariate with the smaller absolute value of correlation
        # to the outcome.
        # Store the indices of the covariates that have to be deleted.
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
        # In case of advanced deletion, delete the covariate with the smaller absolute value of correlation
        # to the outcome, but just if the other covariate contained in the pair has not been deleted yet.
        # Store the indices of the covariates that have to be deleted.
        if ((pos_ordered[i,1] != pos_ordered[i,2]) && !(pos_ordered[i,1] %in% indices_to_delete_with_duplicates)) {
          indices_to_delete_with_duplicates <- append(indices_to_delete_with_duplicates, pos_ordered[i,2])
        }
      }
    }
    # Remove duplicates
    indices_to_delete <- indices_to_delete_with_duplicates[!duplicated(indices_to_delete_with_duplicates)]
    if (length(indices_to_delete)>0) {
      # Return the remaining covariates after deletion
      return(Z[,-indices_to_delete, drop=FALSE])
    } else {
      return(Z)
    }
  } else {
    return(Z)
  }
}

#' Compute interaction terms of covariates
#'
#' @param Z The covariate matrix of dimension nxp
#'
#' @return All interaction terms structured as matrix
#' @export
interaction_terms <- function(Z) {
  # Store the number of covariates
  p <- dim(Z)[2]
  if (p==1) { stop("Input is one-dimensional, no interaction terms computed.") }
  
  # Initialize a matrix of dimension n x (p(p-1)/2) to store all of the p(p-1)/2
  # interaction terms
  ia_terms <- matrix(NA,nrow=dim(Z)[1],ncol=p*(p-1)/2)
  colnames(ia_terms) <- rep("a",p*(p-1)/2)
  column_number <- 1
  # Iterate over all pairs of covariates (without considering order)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      # Calculate and store interaction term
      ia_terms[,column_number] <- Z[,i]*Z[,j]
      colnames(ia_terms)[column_number] <- sprintf("IT %s*%s",colnames(Z)[[i]],colnames(Z)[[j]])
      column_number <- column_number+1
    }
  }
  return(ia_terms)
}

#' Compute cross interaction terms of two sets of covariates
#'
#' @param Z1 The covariate matrix 1 of dimension n x p1
#' @param Z2 The covariate matrix 2 of dimension n x p2
#' @param ident String for naming the resulting columns
#'
#' @return All cross interaction terms structured as a matrix
#' @export
cross_interactions <- function(Z1,Z2,ident="") {
  # Store the numbers of covariates
  p1 <- dim(Z1)[2]
  p2 <- dim(Z2)[2]
  # Store the number of data entries
  n <- dim(Z1)[1]
  
  # Initialize a matrix of dimension n x (p1*p2) to store all of the p1*p2
  # cross interaction terms
  ci_terms <- matrix(NA,nrow=n,ncol=p1*p2)
  colnames(ci_terms) <- rep("a",p1*p2)
  column_number <- 1
  # Iterate over all pairs of covariates (without considering order)
  for(i in 1:p1) {
    for(j in 1:p2) {
      # Calculate and store cross interaction term
      ci_terms[,column_number] <- Z1[,i]*Z2[,j]
      colnames(ci_terms)[column_number] <- sprintf("CI %s: %s * %s",ident,colnames(Z1)[[i]],colnames(Z2)[[j]])
      column_number <- column_number+1
    }
  }
  
  return(ci_terms)
}

#' Calculate trigonometric transformations of the covariates
#'
#' @param Z The covariate matrix of dimension nxp
#' @param order The order up to which the transformations should be calculated
#'
#' @return The trigonometric transformations structured as a matrix
#' @export
fourier_basis <- function(Z, order) {
  # Store the number of covariates
  p <- dim(Z)[2]
  
  # Initialize a matrix of dimension n x (2*p*order) to store all of the
  # 2*p*order trigonometric transformations
  trig_transformations <- matrix(NA,nrow=dim(Z)[1],ncol=p*2*order)
  
  colnames(trig_transformations) <- rep("a",p*2*order)
  col_count <- 1
  # Iterate over all covariates
  for (i in 1:p) {
    # Iterate up to order
    for (j in 1:order) {
      # Calculate and store trigonometric transformations
      trig_transformations[,col_count] <- sin(2*pi*j*Z[,i])
      colnames(trig_transformations)[col_count] <- sprintf("FB %d sin %d",i,j)
      col_count <- col_count+1
      trig_transformations[,col_count] <- cos(2*pi*j*Z[,i])
      colnames(trig_transformations)[col_count] <- sprintf("FB %d cos %d",i,j)
      col_count <- col_count+1
    }
  }
  return(trig_transformations)
}

#' Removes a certain amount of covariates with the highest correlations to other
#' covariates. This function is just used to ensure invertibility in some cases
#' as described in the paragraph "Numerical invertibility" of Section 9.1.
#'
#' @param Z The covariate matrix of dimension nxp
#' @param number The number of covariates to be deleted
#'
#' @return Returns the remaining covariates after deletion
#' @export
remove_covs_with_high_correlation <- function(Z, number) {
  if (isTRUE(ncol(Z)>=1)) {
    covariates_to_delete <- c()
    # Calculate correlation matrix and take the absolute value of each entry
    cor_matrix <- abs(cor(Z))
    # First mark covariates that are not defined (NA) for deletion
    covariates_to_delete <- append(covariates_to_delete,
                                   which(colSums(is.na(cor_matrix)) == nrow(cor_matrix)-1))
    # Fill the according rows and columns of the covariates marked for deletion
    # with zeros. This ensures that they are not taken into consideration again.
    cor_matrix[,covariates_to_delete] <- 0
    cor_matrix[covariates_to_delete,] <- 0
    # Calculate how many covariates still have to be deleted.
    left_to_delete <- number - length(covariates_to_delete)
    if (left_to_delete > 0) {
      # Fill the lower triangle of the correlation matrix with zeros such that
      # we just consider the pairs without taking the order into account
      for (i in 1:nrow(cor_matrix)) {
        for (j in 1:i) {
          cor_matrix[i,j] <- 0
        }
      }
      for (i in 1:left_to_delete) {
        # Identify the entry of the pair with the largest absolute value of correlation
        index <- which(cor_matrix == max(cor_matrix), arr.ind = TRUE)
        # Mark the second covariate contained in that pair for deletion
        covariates_to_delete <- append(covariates_to_delete, index[[1,2]])
        # Fill the respective row and column with zeros so that this covariate
        # is not taken into consideration again
        cor_matrix[index[1,2],] <- 0
        cor_matrix[,index[1,2]] <- 0
      }
    }
    # Return the remaining covariates after deletion
    return(Z[,-covariates_to_delete])
  } else {
    return(Z)
  }
}
