require("mniw")
#' Note cluster counts must be determined before running this procedure
#' Samples a list of matrices denoting the covariance matrices for each cluster
#' 
#' @param A 'n x q' matrix where each row is an observation and the columns are a subset of a larger parameter set.
#' @param mu An list of K 'p x 1' vectors denoting the column means of a cluster.
#' @param N A 'K x 1' vector denoting the counts in each cluster
#' @param K A scalar denoting the number of clusters specified
#' @return An 'K x 1' list of matrices sampled from the inverse-wishart distribution where the kth matrix
#' denotes the posterior covariance matrix of the kth cluster   
#' @details The list of matrices that is returned is  `Sigma_k | A \ \{theta_k\}` which had been sampled from
#' ```
#' `Inverse-Wishart(\sum{\Theta_i - mu_k},N - p - 1)` where N is a vector containing cluster counts
#' , p is a vector of the sizes of theta and 1 is a vector of 1's. Assume that the dimension requirements
#' are met
#' ```

OneStepSampleSigma <- function(theta, mu,N,K){
  n <- length(theta)
  p <- rep(ncol(theta[1],K)) # vector of scalars where each scalar is the size of any theta_i
  theta_mu <- theta - mu # get scale parameter
  # get vector of degrees of freedom
  nu <- N - p - 1
  Psi <- theta_mu %*% theta_mu
  riwish(K,Psi = Psi,nu = nu)
}

