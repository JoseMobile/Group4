require("mniw")

# mu has the same dimension as theta

#' Conditional update for `mu_k` by randomly sampling from the Multivariate Normal model with mean theta and variance proportional to Sigma (see details).
#' 
#' @param theta Estimated random effects of the mixture-normal model. A `n x p` matrix. 
#' @param Sigma Variance matricies of theta. A `p x p x K` array.
#' @param z Estimated classes of the y values. A `n x 1` vector.
#' @param old_mu Current estimated mu. A `K x p` matrix. 
#' @details `K` is the number of classes. `p` is the dimension of one observations of `theta_i`, i.e. the number of features. `n` is the number of observations. If a class `k` does not appear in the z vector, the `k`th row of the return matrix is the `k`th row of currMu. If a class has no observations in the z vector this mu is not updated because having a count of zero causes a division by zero.
#' @return A `K x 1` vector of randomly sampled `mu_k` values from the distribution `mu_k ~ Normal(theta_k, Sigma_k/N_k)`
#' ```
#' where `theta_k` is the mean of all `theta_i`s where `z_i = k` and `N_k` is the number of `z_i = k` (sum of indicator `z_k = k`).
#' 
update_mu <- function(theta, Sigma, z, currMu) {
  # number of classes, K. number of features, q
  SigmaDim <- dim(Sigma)
  K <- SigmaDim[3]  # qxqx*K*
  q <- SigmaDim[1]  # *q*xqxK
  
  # allocating space for class mean theta to be used
  classMeanTheta <- matrix(nrow = K, ncol = q) # Means, theta_k
  
  # if there are no observations of a class k in z, the kth row of Mu won't change
  # after the for loop Mu[k,] is set to currMu[k,]
  # use boolean variable to simplify the process
  replaceMu <- rep(FALSE, K) # initalize to false, changed in for loop if needed
  
  # iterate through the classes
  for (kk in 1:K) {
    zInd <- (z == kk)  # indicator variable if z_i = k
    # have to check how many observations of class in z
    count <- sum(zInd)
    
    # Class mean theta has rows of classMeanTheta filled with mean of class
    if (count > 1) { # more than one observation of class
      classMeanTheta[kk,] <- colMeans(theta[zInd,])
    } else if (count == 1) { # only one observation of class
      classMeanTheta[kk,] <- theta[zInd,]
    } else { # no observations of class
      classMeanTheta[kk,] <- rep(0, q) # fill in theta so that can sample
      replaceMu[kk] <- TRUE # boolean variable used after samplin to replace row
    }
    
    # Calculating scaled Sigma, divided by count of kk in z
    Sigma[,,kk] <- Sigma[,,kk]/count
  }
  
  Mu <- rmNorm(K, classMeanTheta, Sigma) # sample Mu
  
  # replace rows with no observations in that class with last row mean
  Mu[replaceMu] <- currMu[replaceMu] 
  
  return(Mu)
}
