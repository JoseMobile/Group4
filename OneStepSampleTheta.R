require("mniw")

# Let %% denote matrix-vector or matrix-matrix multiplication

#' Samples a matrix, theta_i | A \ {theta_i}, from distribution N(G_i %% (y_i - mu_zi) - mu_zi,G_i %% V_i)
#' where G_i = Sigma_zi %% (V_i - Sigma_zi)^(-1), zi is a randomly sampled index term
#' 
#' @param y An 'n x 1' row vector denoting the features of a single observation.
#' @param V Observed fisher information
#' @param mu An 'n x 1' row vector denoting the prior means.
#' @param Sigma An 'n x n' matrix denoting the prior variances.
#' @return An 'n x 1' matrix that is sampled from mixture normal random effects model.
#' @details The vector that is returned is  `theta_i | A \ \{theta_i\}` which had been sampled from
#' ```
#' `N(G_i(y_i - mu_z_i) - mu_z_i,G_i V_i)``
#' where `G_i = Sigma_z_i(V_i - Sigma_z_i)^(-1)`, and  `z_i` is a randomly sampled index term
#' ```

OneStepSampleTheta <- function(y,V,mu,Sigma){
  theta_i <- rRxNorm(1,y,V,mu, Sigma)
  theta_i
}


# Using slide codes for efficient matrix inversion

#' Solve method for variance matrices.
#'
#' @param V Variance matrix
#' @param x Optional vector or matrix for which to solve system of equations.  If missing calculates inverse matrix.
#' @param ldV Optionally compute log determinant as well.
#' @return Matrix solving system of equations and optionally the log-determinant.
#' @details This function is faster and more stable than `base::solve()` when `V` is known to be positive-definite.
solveV <- function(V, x, ldV = FALSE) {
  C <- chol(V) # cholesky decomposition
  if(missing(x)) x <- diag(nrow(V))
  # solve is O(ncol(C) * ncol(x)) with triangular matrices
  # using backward subsitution
  ans <- backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
  if(ldV) {
    ldV <- 2 * sum(log(diag(C)))
    ans <- list(y = ans, ldV = ldV)
  }
  ans
}

#' Samples group memberships from a multinomial distribution conditioned on the other parameters 
#' 
#' @param theta An 'n x p' matrix, where row i is the parameters of observation i.
#' @param mu An 'K x p' matrix where the row k is expected value of theta if theta is in class k.
#' @param rho An 'K x 1' vector of class membership probabilities 
#' @param Sigma A'p x p x K' array of variance-covariance matrices for each class k.
#' 
#' @return An 'n x 1' vector of updated class memberships.
#' 
#' @details 
update_Z <- function(theta, mu, rho, Sigma){
  # Initialize constants and pre-allocate memory
  K <- length(rho) 
  N <- nrow(theta)
  Kappa <- matrix(NA, nrow=N, ncol=K) 
  
  # Compute inverses and log-determinants for each covariance matrix
  inv_sigma <- apply(Sigma, MARGIN=3, FUN=function(z){ solveV(z, ldV=TRUE) })
  log_det <- sapply(inv_sigma, FUN=function(z){ z$ldV }) 
  log_rho <- log(rho)
  
  # Compute the Kappa matrix 
  for (i in 1:N){
    for (k in 1:K){
      temp <- as.matrix(theta[i, ] - mu[k, ])
      Kappa[i, k] <- log_rho[k] - (t(temp) %*% inv_sigma[[k]]$y %*% temp + log_det[k]) / 2
    }
  }
  
  # Compute the matrix of multinomial probabilties
  Lambda <- exp(Kappa)
  Lambda <- Lambda / rowSums(Lambda) 
  
  # Draw from multinomial N times, each draw of size 1, using different probabilities 
  Z <- apply(Lambda, MARGIN=1, FUN=function(lambda){ rmultinom(1, 1, lambda)} ) # N x K matrix
  # Extract class number using index from matrix of multinomial samples 
  Z <- apply(Z, MARGIN=1, FUN=function(z){ which(z != 0) }) 
  return(Z)
}
