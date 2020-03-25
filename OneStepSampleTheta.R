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
  
#' Samples memberships probabilities from a Dirichlet distribution 
#' 
#' @param theta An 'n x p' matrix, where row i is the parameters of observation i.
#' @param mu An 'n_cat x p' matrix where the row k is expected value of theta if theta is in class k.
#' @param Z An 'n x 1' vector of class memberships for each observation
#' @param Sigma An a list of length 'n_cat', of 'p x p' variance-covariance matrices for each class k.
#' 
#' @return An 'n x 1' vector of updated class membership probabilities.
#' 
#' @details 
update_rho <- function(theta, Z, mu, Sigma){
  # Initialize constants and pre-allocate memory
  K <- nrow(mu) 
  N <- nrow(theta)
  
  # Create alpha vector of class counts
  alpha <- numeric(length=K)
  for (k in 1:K){
    alpha[k] <- sum(Z == k) + 1
  }
  
  # Generate tau's 
  tau <- sapply(alpha, FUN=function(a){ rgamma(1, a, 1) })
  rho <- tau / sum(tau)
  return(rho)
}
