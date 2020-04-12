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




#' Conditional update of theta using the random effects multivariate normal
#' 
#' @param y A `n x p` matrix. Each row of y is an observation on p dimensions of the data
#' @param V A `p x p x n` array of covariance matricies for each observation of y.
#' @param mu A `K x p` matrix. The k-th row of mu is the mean of theta in class k. 
#' @param Sigma A `p x p x K` array of covariance matricies. The k-th matrix is the covariance matrix for theta in class k.
#' @param z A lenght `n` vector for the class observations of y.
#' @details `n` is the number of observations of y. `p` is the dimension of each theta, i.e. the number features of the data. `K` is the number of possible classes of y. 
#' @return A `n x p` matrix. A random sample from the random effects multivariate normal distribution. If the return value is called theta then the details are as follows:
#' theta ~ N(mu, Sigma)
#' y | theta ~ N(theta, V)
#'
#'
update_theta <- function(y, V, mu, Sigma, z) {
  # dimension variables
  K <- nrow(mu) # classes
  n <- nrow(y)  # observations
  p <- ncol(y)  # features
  
  # create new mu and Sigma objects that match 
  # with y over the class observations
  
  # allocate space for new mu and Sigma
  new_mu <- matrix(nrow = n, ncol = p)
  new_Sigma <- array(dim = c(p,p,n))
  
  for (k in 1:K) {
    # iterate through classes and fill in new mu and Sigma with
    # rows of mu and Sigma based on observations from z
    
    ind <- which(z == k) # indices for z that are equal to k
    count <- length(ind) # number of observations of class k
    # fill in all of new_mu with class mean k with index in ind
    new_mu[ind,] <- matrix(mu[k,], nrow = count, ncol = p, byrow = TRUE)
    # fill in indicies of ind of new_Sigma with covariance matricies for k
    new_Sigma[,,ind] <- Sigma[,,k]
  }
  
  # sample using efficient random effects sampler with new mu and Sigma
  theta <- rRxNorm(n, y, V, new_mu, new_Sigma)
  return(theta)
}



#' Conditional update for `mu_k` by randomly sampling from the Multivariate Normal model with mean theta and variance proportional to Sigma (see details).
#' 
#' @param theta Estimated random effects of the mixture-normal model. A `n x p` matrix. 
#' @param old_mu Current estimated mu. A `K x p` matrix. 
#' @param Sigma Variance matricies of theta. A `p x p x K` array.
#' @param z Estimated classes of the y values. A `n x 1` vector.
#' @details `K` is the number of classes. `p` is the dimension of one observations of `theta_i`, i.e. the number of features. `n` is the number of observations. If a class `k` does not appear in the z vector, the `k`th row of the return matrix is the `k`th row of currMu. If a class has no observations in the z vector this mu is not updated because having a count of zero causes a division by zero.
#' @return A `K x 1` vector of randomly sampled `mu_k` values from the distribution `mu_k ~ Normal(theta_k, Sigma_k/N_k)`
#' ```
#' where `theta_k` is the mean of all `theta_i`s where `z_i = k` and `N_k` is the number of `z_i = k` (sum of indicator `z_k = k`).
#' 
update_mu <- function(theta, currMu, Sigma, z) {
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



#' Samples an array of matrices denoting the covariance matrices for each cluster
#' 
#' @param theta A `n x p` matrix where each row is an observation and the columns are a subset of a larger parameter set.
#' @param mu An `K x p` vector denoting the prior cluster means.
#' @param z A length `n` vector denoting the class that each row of theta belongs to.
#' @return A `p x p x K` array of matrices sampled from the inverse-wishart distribution where the kth matrix
#' denotes the posterior covariance matrix of the kth cluster.
#' @details The array of matrices that is returned is  `Sigma_k | A \ \{theta_k\}` which had been sampled from
#' ```
#' `Inverse-Wishart(\sum{\Theta_i - mu_k},N - p - 1)` where N is a vector containing cluster counts
#' , p is a vector of the sizes of theta and 1 is a vector of 1's. Assume that the dimension requirements
#' are met
#' ```
update_Sigma <- function(theta, mu, z) {
  # iterate through classes and make a new mu to match theta
  dim_mu <- dim(mu) # dimensions of mu c(K, p)
  K <- dim_mu[1] # rows are K
  p <- dim_mu[2] # columns are p
  n <- nrow(theta) # rows are n, number of observations
  
  expand_mu <- matrix(nrow = n, ncol = p) # allocate space for expanded mu
  count <- rep(NA, K) # allocate space for count of classes in z
  for (k in 1:K) {
    ind <- which(z == k) # indices of observations of class k
    # calculate count, number of class observations in z
    count[k] <- length(ind)
    
    # fill in expand_mu with each class mean from mu
    expand_mu[ind,] <- matrix(mu[k,], nrow = count[k], ncol = p, byrow = TRUE)
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
    Psi[,,k] <- matrix(0,p,p) # initiate Psi_k to zeros
    for (ii in 1:count[k]) { # iterate through all of the class observations
      ind <- k_inds[ii] # ii-th index of class k
      Psi[,,k] <- Psi[,,k] + res[ind,] %*% t(res[ind,]) # outer product
    }
    # Psi_k is now the sum of outer products of the rows of res
  }
  
  # calculate nu
  nu <- count - p - 1 # a vector of length K
  
  Sigma <- riwish(K, Psi, nu)
  
  return(Sigma)
}



#' Samples cluster memberships probabilities from a Dirichlet distribution 
#' 
#' @param theta An 'n x p' matrix, where row i is the parameters of observation i.
#' @param mu An 'K x p' matrix where the row k is expected value of theta if theta is in class k.
#' @param Sigma A 'p x p x K' array of variance-covariance matrices for each class k.
#' @param Z An 'n x 1' vector of class memberships for each observation
#' @return An 'n x 1' vector of updated class membership probabilities.
update_rho <- function(theta, mu, Sigma, Z){
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



#' Samples group memberships from a multinomial distribution conditioned on the other parameters 
#' 
#' @param theta An 'n x p' matrix, where row i is the parameters of observation i.
#' @param mu An 'K x p' matrix where the row k is expected value of theta if theta is in class k.
#' @param Sigma A'p x p x K' array of variance-covariance matrices for each class k.
#' @param rho An 'K x 1' vector of class membership probabilities 
#' @param give.Lambda A T/F value indicating whether Lambda matrix of probabilities should be returned as well
#' @return An 'n x 1' vector of updated class memberships, or if give.Lambda=TRUE, returns a list that contains the 'n x 1' updated class memberships and the computed 'n x k' Lambda matrix.
update_z <- function(theta, mu, Sigma, rho, give.Lambda=FALSE){
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
  
  # Keep re-sampling until at least one observation in each cluster
  while(TRUE){
    # Draw from multinomial N times, each draw of size 1, using different probabilities 
    Z <- apply(Lambda, MARGIN=1, FUN=function(lambda){ rmultinom(1, 1, lambda)} ) # N x K matrix
    # Extract class number using index from matrix of multinomial samples 
    Z <- apply(Z, MARGIN=1, FUN=function(z){ which(z != 0) }) 

    # Extract class counts 
    counts <- table(Z) 
    if (!(any(counts == 0))){
      break
    }
  }
  
  if (give.Lambda){
    return(list(z=Z, Lambda=Lambda))
  }
  else{
    return(Z)
  }
}
```


```{r}
#' Performs Gibbs Sampling given desired number of clusters and data
#' 
#' @param data The 'N x p' observations for which we want to cluster
#' @param V The 'pxp' observed Fisher information matrix 
#' @param M A scalar value indicating the desired number of MCMC samples to be generated
#' @param K A scalar value indicating the desired number of clusters 
#' @param Theta.out A true or false indicating whether to save all Theta samples
#' @param z.out a true or false indicating whether to save all z samples 
#' @return A list of all the sampled 
gibbsSampler <- function(data, V, M=1000, K, Theta.out=FALSE, z.out=FALSE){
  N <- nrow(data)
  p <- ncol(data)
  
  # Generate Initial values 
  old_theta <- data
  old_mu <- array(rnorm(K*p), nrow=K, ncol=p)
  old_Sigma <- array(rnorm(p*p*K), dim=c(p, p, K))
  for (k in 1:K){
    sigma_0[, , k] <- t(sigma_0[, , k]) %*% sigma_0[, , k]
  }
  old_rho <- rep(1/K, times=K)
  old_z <- sample(1:K, N, replace=TRUE)

  # Allocate space to store all generated samples
  rho <- array(NA, dim=c(K, M))
  mu <- array(NA, dim=c(K, p, M))
  Sigma <- array(NA, dim=c(p, p, K, M))
  Theta <- ifelse(Theta.out, array(NA, dim=c(p, N, M)), NULL)
  z <- ifelse(z.out, array(NA, dim=c(N, M)), NULL)
  
  total_Lambda <- array(0, dim=c(N, K))
  
  # The dimensions of the arguments are: 
  #   theta: N x p      mu: K x p     Sigma: p x p x K    rho: K x 1    z: N x 1
  for (m in 1:M){
    # Call the update functions
    #  Check the dimensions of the return values match the other functions 
    
    # Update theta, a 'N x p' matrix of random effects of y
    new_theta <- update_theta(data, V, old_mu, old_Sigma, old_z)
    
    # Update mu, a 'K x p' matrix of class averages of theta
    mu[, , m] <- update_mu(old_theta, old_mu, old_Sigma, old_z)
    
    # Update Sigma, a 'p x p x K' array of covariance matricies
    Sigma[, , , m] <- update_Sigma(old_theta, old_mu, old_z)
    
    # Update rho, a K x 1 vector
    rho[, m] <- update_rho(old_theta, old_z, old_mu, old_Sigma)
    
    # Update z, an N x 1 vector
    new_z <- update_z(old_theta, old_mu, old_rho, old_Sigma, give.Lambda=TRUE)
    
    # Store theta and z values if desired
    if (Theta.out){
      Theta[, , m] <- new_theta
    }
    if (z.out){
      z[, m] <- new_z
    }
    
    # Update all values
    old_theta <- new_theta
    old_mu <- mu[, , m]
    old_Sigma <- Sigma[, , , m]
    old_rho <- rho[, m]
    old_z <- new_z 
    
    total_Lambda <- total_Lambda + new_z$Lambda
  }
  
  total_Lambda <- total_Lambda / M

  return(list(theta=Theta, mu=mu, Sigma=Sigma, rho=rho, z=z, Lambda=total_Lambda))
}
