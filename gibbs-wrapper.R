source('nnm-functions.R')



#' Conditional update of theta using the random effects multivariate normal
#' 
#' 
#' @param y A `n x p` matrix. Each row of y is an observation on p dimensions of the data
#' @param V A `p x p x n` array of covariance matricies for each observation of y.
#' @param mu A `K x p` matrix. The k-th row of mu is the mean of theta in class k. 
#' @param Sigma A `p x p x K` array of covariance matricies. The k-th matrix is the covariance matrix for theta in class k.
#' @param z A lenght `n` vector for the class observations of y.
#' @details `n` is the number of observations of y. `p` is the dimension of each theta, i.e. the number features of the data. `K` is the number of possible classes of y. 
#' 
#' @return A `n x p` matrix. A random sample from the random effects multivariate normal distribution. If the return value is called theta then the details are as follows:
#' theta ~ N(mu, Sigma)
#' y | theta ~ N(theta, V)
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
  count <- numeric(length = K)

  
  for (k in 1:K) {
    # iterate through classes and fill in new mu and Sigma with
    # rows of mu and Sigma based on observations from z
    
    ind <- which(z == k) # indices for z that are equal to k
    count[k] <- length(ind) # number of observations of class k
    # fill in all of new_mu with class mean k with index in ind
    # only execute if there are class observations to avoid warnings
    if (count[k] > 0) { 
      new_mu[ind,] <- matrix(mu[k,], nrow = count[k], ncol = p, byrow = TRUE)
      # fill in indicies of ind of new_Sigma with covariance matricies for k
      new_Sigma[,,ind] <- Sigma[,,k]
    }
  }
  
  # Checking the dimensions of y and new_Sigma
  dim_nSig <- dim(new_Sigma) # dimensions of new sigma
  dim_y <- dim(y)            # dimensions of y
  
  # if statement for if anything in the function looks weird
  if (dim_y[1] != dim_nSig[3] || # n in y is not the same as n in Sigma
      dim_y[2] != dim_nSig[1] || # p in y is not the same as p in Sigma
      dim_y[2] != dim_nSig[2] || # p in y is not the same as the other p in Sigma
      anyNA(new_Sigma)) 
  { 
    cat("dimensions of y :        ", dim_y, "\n")
    cat("dimensions of new_Sigma :", dim_nSig, "\n")
    
    cat("any NAs in y :", anyNA(y), "\n")
    cat("any NAs in new Sigma :", anyNA(new_Sigma), "\n")
    cat("count of classes :", count, "\n")
    cat("any NAs in new mu :", anyNA(new_mu), "\n")
    
    cat("any NULLs in y :", is.null(y), "\n")
    cat("any NULLs in new Sigma :", is.null(new_Sigma), "\n")
    cat("any NULLs in new mu :", is.null(new_mu), "\n")
    
    print(Sigma)
    print("the above is theSigma array")
  }
  
  # sample using efficient random effects sampler with new mu and Sigma
  theta <- rRxNorm(n, y, V, new_mu, new_Sigma)
  return(theta)
}



#' Conditional update for `mu_k` by randomly sampling from the Multivariate Normal model with mean theta and variance proportional to Sigma (see details).
#' 
#' 
#' @param theta Estimated random effects of the mixture-normal model. A `n x p` matrix. 
#' @param currMu Current estimated mu. A `K x p` matrix. 
#' @param Sigma Variance matricies of theta. A `p x p x K` array.
#' @param z Estimated classes of the y values. A `n x 1` vector.
#' 
#' @details `K` is the number of classes. `p` is the dimension of one observations of `theta_i`, i.e. the number of features. `n` is the number of observations. If a class `k` does not appear in the z vector, the `k`th row of the return matrix is the `k`th row of currMu. If a class has no observations in the z vector this mu is not updated because having a count of zero causes a division by zero.
#' @return A `K x 1` vector of randomly sampled `mu_k` values from the distribution `mu_k ~ Normal(theta_k, Sigma_k/N_k)`
#' ```
#' where `theta_k` is the mean of all `theta_i`s where `z_i = k` and `N_k` is the number of `z_i = k` (sum of indicator `z_k = k`).
#' 
update_mu <- function(theta, old_mu, Sigma, z) {
  # number of classes, K. number of features, q
  SigmaDim <- dim(Sigma)
  K <- SigmaDim[3]  # q x q x K
  q <- SigmaDim[1]  # q x q x K
  
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
      
      # Calculating scaled Sigma, divided by count of kk in z
      Sigma[,,kk] <- Sigma[,,kk]/count
    } else if (count == 1) { # only one observation of class
      classMeanTheta[kk,] <- theta[zInd,]
      
      # Calculating scaled Sigma, divided by count of kk in z
      Sigma[,,kk] <- Sigma[,,kk]/count
    } else { # no observations of class
      classMeanTheta[kk,] <- rep(0, q) # fill in theta so that can sample
      replaceMu[kk] <- TRUE # boolean variable used after samplin to replace row
    }
    

  }
  
  Mu <- rmNorm(K, classMeanTheta, Sigma) # sample Mu
  
  # replace rows with no observations in that class with last row mean
  Mu[replaceMu,] <- old_mu[replaceMu,] 
  
  return(Mu)
}


#' Samples an array of matrices denoting the covariance matrices for each cluster
#' 
#' 
#' @param y A `n x p` matrix. Each row of y is an observation on p dimensions of the data.
#' @param theta A `n x p` matrix where each row is an observation and the columns are a subset of a larger parameter set.
#' @param mu An `K x p` vector denoting the prior cluster means.
#' @param z A length `n` vector denoting the class that each row of theta belongs to.
#' @param Omega The covariance matrixes of our chosen prior. 
#' @param vK The virtual counts corresponding to mu of chosen prior
#' 
#' @return A `p x p x K` array of matrices sampled from the inverse-wishart distribution where the kth matrix
#' denotes the posterior covariance matrix of the kth cluster.
#' @details The array of matrices that is returned is  `Sigma_k | A \ \{theta_k\}` which had been sampled from
#' ```
#' `Inverse-Wishart(\sum{\Theta_i - mu_k},N + vK)` where N is a vector containing cluster counts
#' , p is a vector of the sizes of theta and 1 is a vector of 1's. Assume that the dimension requirements
#' are met
#' ```
#' 
update_Sigma <- function(y, theta, mu, old_Sigma, z, Omega, vK) {
  # iterate through classes and make a new mu to match theta
  dim_mu <- dim(mu) # dimensions of mu c(K, p)
  K <- dim_mu[1] # rows are K
  p <- dim_mu[2] # columns are p
  n <- nrow(theta) # rows are n, number of observations
  
  # allocate space before 
  replaceSigma <- rep(FALSE, K) # indicator for if there are empty clusters
  # initialized at false, will change if empty cluster if found
  count <- rep(NA, K) # allocate space for count of classes in z
  Psi <- array(dim = c(p,p,K)) # allocate space for Psi matrix
  
  for (k in 1:K) {
    k_inds <- which(z == k) # indices of observations of class k
    # calculate count, number of class observations in z
    count[k] <- length(k_inds)
    
    # initiate Psi_k to Omega_k
    Psi[,,k] <- Omega[,,k] 
    
    # check how many observations there are in class k
    if (count[k] == 0) {
      replaceSigma[k] <- TRUE # found empty cluster
      print("cluster of 0")
    } else { # cluster has at least one observation
      for (ii in 1:count[k]) { # iterate through all of the class observations
        ind <- k_inds[ii] # ii-th observation of class k
        
        # calculate residual of theta and mu
        res <- theta[ind,] - mu[k,]
        Psi[,,k] <- Psi[,,k] + tcrossprod(res) # outer product
        # tcrossprod(x) is x %*% t(x)
      }
      # Psi_k is now the sum of outer products of the rows of res
    }
  }
  
  # calculate nu
  nu <- count + vK # a vector of length K
  Sigma <- riwish(K, Psi, nu)
  
  # replace empty cluster class Sigmas with the previous Sigma
  Sigma[,,replaceSigma] <- old_Sigma[,,replaceSigma]
  
  return(Sigma)
}

#' Samples cluster memberships probabilities from a Dirichlet distribution 
#' 
#' 
#' @param theta An 'n x p' matrix, where row i is the parameters of observation i.
#' @param mu An 'K x p' matrix where the row k is expected value of theta if theta is in class k.
#' @param Sigma A 'p x p x K' array of variance-covariance matrices for each class k.
#' @param Z An 'n x 1' vector of class memberships for each observation
#' 
#' @return An 'n x 1' vector of updated class membership probabilities.
#' 
update_rho <- function(theta, mu, Sigma, Z){
  # Initialize constants and pre-allocate memory
  K <- nrow(mu) 
  N <- nrow(theta)
  
  # Create alpha vector of class counts
  alpha <- numeric(length=K)
  for (k in 1:K){
    alpha[k] <- sum(Z == k) + 1
  }
  
  rho<-rdirichlet(1, alpha)
  return(rho)
}


#' Samples group memberships from a multinomial distribution conditioned on the other parameters 
#' 
#' 
#' @param theta An 'n x p' matrix, where row i is the parameters of observation i.
#' @param mu An 'K x p' matrix where the row k is expected value of theta if theta is in class k.
#' @param Sigma A'p x p x K' array of variance-covariance matrices for each class k.
#' @param rho An 'K x 1' vector of class membership probabilities 
#' @param give.Lambda A T/F value indicating whether Lambda matrix of probabilities should be returned as well
#' 
#' @return An 'n x 1' vector of updated class memberships, or if give.Lambda=TRUE, returns a list that contains the 'n x 1' updated class memberships and the computed 'n x k' Lambda matrix.
#' 
update_z <- function(theta, mu, Sigma, rho){
  # Initialize constants and pre-allocate memory
  K <- length(rho) 
  N <- nrow(theta)
  Kappa <- matrix(NA, nrow=N, ncol=K) 
  
  # Compute inverses and log-determinants for each covariance matrix
   inv_sigma <- vector("list", length=K)
   for (k in 1:K){
     inv_sigma[[k]] <- solveV(Sigma[, , k], ldV=TRUE)
   }
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

  Z <- rcategorical(t(Lambda))
  
  return(list(z=Z, Lambda=Lambda))
}

###########################################################################################3



#' Performs Gibbs Sampling given desired number of clusters and data
#' 
#' @param data The `N x p`` observations for which we want to cluster
#' @param V The `p x p x N` observed Fisher information matrix 
#' @param prior A list containing the parameters for the prior selected 
#' @param initParamVals A list containing the initial values of the parameters of the model
#' @param K A scalar value indicating the desired number of clusters 
#' @param burnin_period A scalar value indicating the initial number of samples to discard
#' @param numIter scalar value indicating the desired number of MCMC samples to be generated
#' @param mu.fixed A true or false indicating whether to sample new mu values 
#' @param Sigma.fixed  A true or false indicating whether to sample new Sigma values 
#' @param theta.fixed A true or false indicating whether to sample new theta values 
#' @param rho.fixed A true or false indicating whether to sample new rho values 
#' @param z.fixed A true or false indicating whether to sample new z values 
#' @param Theta.out A true or false indicating whether to save all Theta samples
#' @param z.out a true or false indicating whether to save all z samples 
#' @return A list of all the sampled 
gibbsSampler <- function(data, V, prior, initParamVals, K, burnin_period, numIter, mu.fixed = FALSE, 
                         Sigma.fixed =FALSE, theta.fixed= FALSE, rho.fixed = FALSE, z.fixed = FALSE, Theta.out=FALSE, z.out=FALSE){
  N <- nrow(data)
  p <- ncol(data)
    
  #Getting the Prior Parameters
  if (missing(prior)){
    Vk <- rep(p+2,K)
    Omega<- array(rep(0, p*p*K), dim=c(p,p,K))
    
    for(i in 1:K)
      Omega[,,i]<- var(data) *(N-1)/(N-p)
  } else {
    if ( (! "Vk" %in% names(prior) ) || (!"Omega" %in% names(prior)) )
      stop("prior needs Vk and Omega params!")
    
    Vk<- prior$Vk 
    Omega<- prior$Omega
  }
  
  #Getting the Initial Value of Model Parameters 
  if (missing(initParamVals)){ # no initial values given
    # INITALIZING PARAMETERS
    
    # Compute overall data means and variances
    data_mu <- colMeans(data)
    data_var <- var(data) * (N-1)/(N-p) # variance of data
    
    # initalizing parameters using kmeans++ algorithm
    old_z <- init_z(data, K) # use kmeans++ algorith to initialize z
    
    # set rho to proportion of each observation
    count <- as.numeric(table(old_z))
    
    # table gives count of each observation of class
    # then dividing by N gives the proportion
    old_rho <- count/N
    
    # set theta to observed data
    old_theta <- data
    
    # allocate space for mu and Sigma
    old_mu <- matrix(nrow = K, ncol = p)
    old_Sigma <- array(dim = c(p,p,K))
    # set mu to sample mean of data in that class
    for (k in 1:K) {
      ind <- which(old_z == k) # indices of data in class k
      if (count[k] > p) { # enough data to sample both mu and Sigma
        old_mu[k,] <- colMeans(data[ind,]) # sample mu
        # sample Sigma, var() scaled by n-1, so rescale to n-p
        old_Sigma[,,k] <- var(data[ind,]) * (count[k]-1)/(count[k]-p)
      } else {
        # not enough data for Sigma, so use data variance
        old_Sigma[,,k] <- data_var
        # sampling for mu
        if (count[k] > 1) { # enough data for sampling mu
          old_mu[k,] <- colMeans(data[ind,])
        } else if (count[k] == 1) { # only one so set to mu
          old_mu[k,] <- data[ind,]
        } else { # count is zero, so use data mean
          old_mu[k,] <- data_mu
        } # end second if
      } # end first if
    } # end for
    
  } else { # given inital parameters
    if ( (! "theta" %in% names(initParamVals) ) || (!"mu" %in% names(initParamVals)) ||  (!"Sigma" %in% names(initParamVals))
         ||  (!"rho" %in% names(initParamVals))||  (!"z" %in% names(initParamVals)))
      stop("initParamVals needs theta, mu, Sigma, rho and z params!")
    
    old_theta<- initParamVals$theta
    old_mu<- initParamVals$mu
    old_Sigma<- initParamVals$Sigma
    old_rho<- initParamVals$rho
    old_z<- initParamVals$z
  }

  # Allocate space to store all generated samples
  rho <- array(NA, dim=c(K, numIter))
  mu <- array(NA, dim=c(K, p, numIter))
  Sigma <- array(NA, dim=c(p, p, K, numIter))
  if (Theta.out)
    Theta <- array(NA, dim=c(p, N, numIter))
  else 
    Theta<- NULL 
  
  if (z.out)
    z <- array(NA, dim=c(N, numIter))
  else 
    z<- NULL
  
  total_Lambda <- array(0, dim=c(N, K))
  
  # The dimensions of the arguments are: 
  #   theta: N x p      mu: K x p     Sigma: p x p x K    rho: K x 1    z: N x 1
  # Parameters orders should all be: theta, mu, Sigma, rho, z
  #   If not check the function calls and standardize them
  for (m in 1:(burnin_period +numIter)){
    # Call the update functions
    
    #Keep Track of Current Iteration Number
    if (mod(m, 200) ==0 ){
      cat("Iteration #", m, "\n")
    }
    
    if (!theta.fixed){
    # Update theta, a 'N x p' matrix of random effects of y
      new_theta <- update_theta(data, V, old_mu, old_Sigma, old_z)
    }
    if (!mu.fixed){
    # Update mu, a 'K x p' matrix of class averages of theta
      new_mu <- update_mu(new_theta, old_mu, old_Sigma, old_z)
    }
    if(!Sigma.fixed){
    # Update Sigma, a 'p x p x K' array of covariance matricies
      new_Sigma <- update_Sigma(data, new_theta, new_mu, old_Sigma, old_z, Omega, Vk)
    }
    if(!rho.fixed){
    # Update rho, a K x 1 vector
      new_rho <- update_rho(new_theta, new_mu, new_Sigma, old_z)
    }
    if (!z.fixed){
    # Update z, an N x 1 vector
      new_z <- update_z(new_theta, new_mu, new_Sigma, new_rho)
    }
    
    # Store theta and z values if desired
    if (m > burnin_period){
      curIdx<- m - burnin_period
      if (Theta.out){
        Theta[, , curIdx] <- new_theta
      }
      if (z.out){
        z[, curIdx] <- new_z$z
      }
      
      rho[ ,curIdx]<- new_rho
      mu[ , ,curIdx]<- new_mu
      Sigma[ , , ,curIdx]<-  new_Sigma
      total_Lambda <- total_Lambda + new_z$Lambda
    }
    
    # Update all values
    old_theta <- new_theta
    old_mu <- new_mu
    old_Sigma <- new_Sigma
    old_rho <- new_rho
    old_z <- new_z$z
  }
  
  total_Lambda <- total_Lambda / numIter
  
  return(list(theta=Theta, mu=mu, Sigma=Sigma, rho=rho, z=z, Lambda=total_Lambda))
}
