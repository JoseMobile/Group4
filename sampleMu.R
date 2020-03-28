require("mniw")

# mu has the same dimension as theta

#' Conditional update for `mu_k` by randomly sampling from the Multivariate Normal model with mean theta and variance proportional to Sigma (see details).
#' 
#' @param theta Estimated random effects of the mixture-normal model. A `n x q` matrix. 
#' @param Sigma Variance matricies of theta. A `q x q x K` array.
#' @param z Estimated classes of the y values. A `n x 1` vector.
#' @details `K` is the number of classes. `q` is the dimension of one observations of `theta_i`, i.e. the number of features. `n` is the number of observations.
#' @return A `K x 1` vector of randomly sampled `mu_k` values from the distribution `mu_k ~ Normal(theta_k, Sigma_k/N_k)`
#' ```
#' where `theta_k` is the mean of all `theta_i`s where `z_i = k` and `N_k` is the number of `z_i = k` (sum of indicator `z_k = k`).
#' 
sampleMu <- function(theta, Sigma, z) {
  # number of classes, K. number of features, q
  SigmaDim <- dim(Sigma)
  K <- SigmaDim[3]  # qxqx*K*
  q <- SigmaDim[1]  # *q*xqxK
  
  # allocating space for class mean theta to be used
  classMeanTheta <- matrix(nrow = K, ncol = q) # Means, theta_k
  
  # iterate through the classes
  for (kk in 1:K) {
    zInd <- ifelse(z == kk, yes = 1, no = 0)  # indicator variable if z_i = k
    
    # Class mean theta has rows of classMeanTheta filled with mean of class
    classMeanTheta[kk,] <- colMeans(theta[zInd,])
    
    # Calculating scaled Sigma, divided by count of kk in z
    Sigma[,,kk] <- sum(z == kk)
  }
  
  Mu <- rmNorm(K, classMeanTheta, Sigma)
  
  return(Mu)
}
