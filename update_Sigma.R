require("mniw")

#' Samples an array of matrices denoting the covariance matrices for each cluster
#' 
#' @param y A `n x p` matrix. Each row of y is an observation on p dimensions of the data.
#' @param theta A `n x p` matrix where each row is an observation and the columns are a subset of a larger parameter set.
#' @param mu An `K x p` vector denoting the prior cluster means.
#' @param z A length `n` vector denoting the class that each row of theta belongs to.
#' @param old_Sigma The current observation of the Sigma. A `p x p x K` array of covariance matrices.
#' @return A `p x p x K` array of matrices sampled from the inverse-wishart distribution where the kth matrix
#' denotes the posterior covariance matrix of the kth cluster.
#' @details The array of matrices that is returned is  `Sigma_k | A \ \{theta_k\}` which had been sampled from
#' ```
#' `Inverse-Wishart(\sum{\Theta_i - mu_k},N - p - 1)` where N is a vector containing cluster counts
#' , p is a vector of the sizes of theta and 1 is a vector of 1's. Assume that the dimension requirements
#' are met
#' ```
update_Sigma <- function(y, theta, mu, z) {
  # iterate through classes and make a new mu to match theta
  dim_mu <- dim(mu) # dimensions of mu c(K, p)
  K <- dim_mu[1] # rows are K
  p <- dim_mu[2] # columns are p
  n <- nrow(theta) # rows are n, number of observations
  
  expand_mu <- matrix(nrow = n, ncol = p) # allocate space for expanded mu
  count <- rep(NA, K) # allocate space for count of classes in z
  Omega <- array(dim = c(p,p,K)) # allocate space for omega
  for (k in 1:K) {
    ind <- which(z == k) # indices of observations of class k
    # calculate count, number of class observations in z
    count[k] <- length(ind)
    
    # fill in expand_mu with each class mean from mu
    expand_mu[ind,] <- matrix(mu[k,], nrow = count[k], ncol = p, byrow = TRUE)
    
    # fill omega array with variances of y
    Omega[,,k] <- var(y) # variance matrix of y
  }
  # now the rows of theta and mu match so that the i-th row of mu is the mean
  # for the i-th row of theta
  
  # residuals from theta and mus
  res <- theta - expand_mu 
  
  # calculate Psi matrix
  Psi <- array(dim = c(p,p,K)) # allocate space for Psi matrix
  # calculate Psi by outer product of each row of res an entry to the array
  for (k in 1:K) {
    k_inds <- which(z == k) # indicies of observations of class k
    Psi[,,k] <- Omega[,,k] # initiate Psi_k to Omega_k
    for (ii in 1:count[k]) { # iterate through all of the class observations
      ind <- k_inds[ii] # ii-th index of class k
      Psi[,,k] <- Psi[,,k] + tcrossprod(res[ind,]) # outer product
      # tcrossprod(x) is x %*% t(x)
    }
    # Psi_k is now the sum of outer products of the rows of res
  }
  
  # calculate nu
  nu <- count + p + 2 # a vector of length K
  
  Sigma <- riwish(K, Psi, nu)
  
  return(Sigma)
}

# correct dimension numbers to quickly check code
#K <- 7
#n <- 1000
#p <- 10
#theta <- matrix(rnorm(n*p), nrow = n)
#mu <- matrix(rnorm(K*p), nrow = K)
#z <- ceiling(K*runif(n))
#y <- matrix(rnorm(n*p), nrow = n)
#update_Sigma(y, theta, mu, z)
