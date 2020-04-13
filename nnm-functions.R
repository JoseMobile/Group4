require(mniw)
require(numbers)

#' Complete data loglikelihood for the normal-normal-mixture (NNM) model.
#'
#'
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param rho A length-`K` probability vector of group membership.
#' @param y An `N x p` matrix of observations, each of which is a column.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#' @param theta A `p x N` matrix of observation means.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param log bool value indicating whether to give log value of distribution 
#' 
#' @return scalar value of full data loglikelihood 
nnm_loglik <- function(mu, Sigma, rho, y, V, theta, z, log=TRUE) {
  K <- nrow(mu) # number of groups
  N <- nrow(y) #number of data points 
  ll <- 0
  for(kk in 1:K) {
    idk <- z == kk # T/F membership in group k
    Nk <- sum(idk) # number of observations in group k

    if(Nk > 0) {
        # group membership contribution
        ll <- ll + Nk * log(rho[kk])
        # mean contribution
        # here we repeat mu[kk,] and Sigma[,,kk] Nk times to correspond to the total number of group membership count  
        ll <- ll + sum(dmNorm(t(theta[,idk]), mu = matrix(rep(mu[kk,], Nk), nrow=Nk, ncol=p, byrow=TRUE), Sigma = array(rep(Sigma[,,kk], Nk), dim=c(p,p,Nk)),  log = TRUE ))
    }
  }
  # observation contribution
  ll<- ll + sum(dmNorm(y, mu = t(theta), Sigma = V, log = TRUE))
  
  if (!log) # check if non log value desired  
    ll<- exp(ll)
  
  ll
}


#' Complete data posterior for the normal-normal-mixture (NNM) model.
#'
#'
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param rho A length-`K` probability vector of group membership.
#' @param y An `N x p` matrix of observations, each of which is a column.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#' @param theta A `p x N` matrix of observation means.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param vk A length- `K` vector of integers acting as virtual counts from prior 
#' @param Omega A `p x p x K` matrix of variances corresponding to the variance of the prior chosen.
#' @param log bool value indicating whether to give log value of distribution 
#' 
#' @return scalar value of full data log posterior 
nnm_post <- function(mu, Sigma, rho, y, V, theta, z, vk, Omega, log=TRUE) {
  K<- nrow(mu) # number of groups 
  
  ll<- nnm_loglik(mu, Sigma, rho, y, V, theta, z) #complete data loglikelihood 
  
  for(i in 1:K){
    ll<- ll - 0.5* ( (vk[i] + p+1)*log(det(Sigma[,,i])) + sum(diag(solveV(Sigma[,,i]) %*% Omega[,,i])) )
  } # complete data log posterior 

  if (!log) # check if non log value desired  
    ll <- exp(ll) 

  ll 
}

#' Normalized Conditional Log-PDF of Mu 
#' 
#' 
#' @param theta A `p x N` matrix of observation means.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param k An integer value beween 1 and `K` to indicate which mu to condition on 
#' @param log bool value indicating whether to give log value of distribution 
#' 
#' @return scalar value of log P(mu_k | A\{mu_k}) 
mu_k.post<-function(theta, mu, Sigma, z, k, log=TRUE){
  #extract mu and sigma of the kth group
  mu_k <- mu[k,] 
  Sigma_k <- Sigma[,,k]
  
  #get number of group participants
  tidx <-  tidx<-(k == z)
  Nk<- sum(tidx)
  
  #if Nk == 1 then no need for rowMeans
  theta_avg_k<- theta[,tidx]  
  
  #Nk > 1, avg group members values across each dimension
  if (Nk > 1){
    theta_avg_k<-rowMeans(theta[,tidx])
  }
  
  #ensure vector data structure is used  
  theta_avg_k<-as.vector(theta_avg_k)
  
  #get log pdf of a multivariate normal distribution 
  ll<-dmNorm(mu_k, theta_avg_k, Sigma_k/Nk, log=TRUE)
  
  if(!log) # check if non log value desired  
    ll<-exp(ll)
  
  ll
}


#' Extract Mean and Variance Marginals of parameter Mu_k assuming stationary theta and Sigma, and z
#' 
#' 
#' @param theta A `p x N` matrix of observation means.
#' @param marginal_idx a scalar value specifying which demension of mu to condition on 
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param k An integer value beween 1 and `K` to ondicate which mu to condition on 
#' 
#' @return list containing both marginal mean and sd of Mu dat  
mu_k_marginals.param<-function(theta, marginal_idx, Sigma, z, k){
  
  Sigma_k <- Sigma[,,k]
  
  #get group size or group membership count 
  tidx <-  tidx<-(k == z)
  Nk<- sum(tidx)
  
  #if Nk == 1 then theta_avg is equal to only member 
  theta_avg_k<- theta[,tidx]  
  
  #Nk > 1, avg group members values across each dimension
  if (Nk > 1){
    theta_avg_k<-rowMeans(theta[,tidx])
  }
  
  Sigma_k <- Sigma_k/Nk
  
  #extract desired marginal params as specified by marginal_idx 
  params<- list(mean= theta_avg_k[marginal_idx], sd=sqrt(Sigma_k[marginal_idx, marginal_idx]))
  
  params
}


#' Normalized Conditional Log-PDF of Theta_i
#' 
#' 
#' @param y An `N x p` vector of observations, each of which is a column.
#' @param theta A `p x N` vector of observation means.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param Z A length- `N` vector of integers between 1 and `K`
#' @param rand_idx scalar value netween `1` and `N` to indicate which theta to condition on 
#' @param log bool value indicating whether to give log value of distribution 
#' 
#' @return scalar value of log P(theta_i| A\{theta_i}) 
theta_i.post<-function(y, theta, V, mu, Sigma, Z, rand_idx, log=TRUE){
  
  #Extract params relevant to the theta conditoning on 
  y_i <- y[rand_idx,]
  theta_i<- theta[,rand_idx]
  V_i<- V[,, rand_idx]
  z_i <- Z[rand_idx]
  
  #get Sigma and mu of group theta_i belongs to (according to z_i) 
  Sigma_zi <- Sigma[,,z_i]
  mu_zi <- mu[z_i,]
  
  #Precompute the new Mu and Variance of the resulting distribution 
  G_i <- Sigma_zi %*% solveV(V_i + Sigma_zi)
  new_mu <- G_i  %*% (y_i - mu_zi) + mu_zi
  new_V <- G_i  %*% V_i 
  
  #ensure that new_mu is a vector 
  new_mu<-as.vector(new_mu)
  
  #get log pdf of a multivariate normal distribution
  ll<- dmNorm(theta_i, new_mu, new_V, log=TRUE)
  
  if(!log)  # check if non log value desired  
    ll<-exp(ll)
  
  ll
}


#' Normalized Conditional Log-PDF of Sigma_k
#' 
#'
#' @param theta A `p x N` matrix of observation means.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param Z A length- `N` vector of integers between 1 and `K`
#' @param k An integer value beween 1 and `K` to indicate which Sigma to condition on 
#' @param vk A length- `K` vector of integers acting as virtual counts from prior 
#' @param Omega  A `p x p x K` matrix of variances corresponding to the variance of the prior chosen.
#' @param log bool value indicating whether to give log value of distribution 
#' 
#' @return scalar value of log P(Sigma_k| A\{Sigma_k}) 
sigma_k.post<-function(theta, mu, Sigma, Z, k, vk, Omega, log=TRUE){
  
  #extract relevant params corresponding to the kth group 
  Sigma_k <- Sigma[,,k] #extract Sigma corresponding to the kth group 
  mu_k <- mu[k, ]
  vkk<- vk[k]
  Omegak <- Omega[,,k]
  
  Sigma_inv<- solveV(Sigma_k) #get the inverse of Sigma_k 
  
  #get total membership counts of kth group 
  tidx<- (Z == k)
  Nk<- sum(tidx)
  
  theta_less_mu <- theta[,tidx] - mu_k
  
  #get log pdf of an inverse-wishart distribution
  ll <- diwish(Sigma_k, Omegak+(theta_less_mu %*% t(theta_less_mu)) ,  Nk + vkk, log = TRUE)
  
  if (!log) # check if non log value desired
    ll<-exp(ll)
  
  ll
}


#' Normalized Conditional Log-PDF of Rho
#' 
#'
#' @param Z A length- `N` vectir of integers between 1 and `K`
#' @param rho A length-`K` probability vector of group membership.
#' @param log bool value indicating whether to give log value of distribution 
#' 
#' @return scalar value of log P(rho| A\{rho}) 
rho.post<-function(Z, rho, log=TRUE){
  #Get total number of groups 
  K<- length(rho)
  
  #Initialize group counts to 0 
  Nk<- rep(0, K)
  
  #Precompute log value of rho elements
  log_rho <- log(rho)
  
  #get group counts for each group 
  for(j in 1:K){
    tidx<- (j == Z)
    Nk[j] <- sum(tidx)
  }
  
  #get log pdf of an Dirichlet distribution
  ll<-ddirichlet(rho, Nk + 1, log=TRUE)
  
  if (!log)  # check if non log value desired
    ll<-exp(ll)
  
  ll
}


#' Normalized Conditional Log-PDF of Z_i 
#' 
#' 
#' @param Z A length- `N` vectir of integers between 1 and `K`
#' @param rho A length-`K` probability vector of group membership.
#' @param theta A `p x N` vector of observation means.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param rand_idx A scalar value bewteen `1` and `N` to indicate which z to condtion on. 
#' @param log bool value indicating whether to give log value of distribution 
#' 
#' @return scalar value of log P(z_i = K | A\{z_i}) 
z.post<-function(Z, rho, theta, mu, Sigma, rand_idx, log=TRUE){
  
  #initialize log posterior value to 0
  ll<-0 
  
  #get group membership theta_i 
  z_i <- Z[rand_idx]
  theta_i <- theta[, rand_idx]
  
  #add rho contribution 
  ll<- ll + log(rho[z_i])
  
  #get mu and Sigma of group theta_i belongs 
  Sigma_k <- Sigma[,, z_i]
  mu_k <- mu[z_i,]
  
  #get log pdf of an multivariate normal distribution
  ll<- ll + dmNorm(theta_i, mu= mu_k, Sigma= Sigma_k, log=TRUE)
  
  if (!log) # check if non log value desired
    ll<-exp(ll)
  
  ll
}


#' Sample from a categorical distribution.
#'
#'
#' Performs one draw from a categorical distribution (i.e., a multinomial distribution with size `n = 1`) for each of multiple probability vectors.
#'
#' @param prob An `n_cat x n_prob` matrix of probability vectors, each of which is a column of `prob`.  The entries of `prob` must be nonnegative, but are internally normalized such that each column sums to one.
#' @return A vector of length `n_prob` of integers between 1 and `n_cat` indicating the category draw for each column of `prob`.
rcategorical <- function(prob) {
  if(any(prob < 0)) stop("prob must contain only nonnegative entries.")
  cprob <- apply(prob, 2, cumsum)
  u <- runif(ncol(prob)) * cprob[nrow(prob),]
  
  apply(sweep(cprob, 2, u) >= 0, 2, which.max)
}


#' Density of the Dirichlet distribution.
#'
#'
#' @param x Observation vector of nonnegative entries which sum to one.
#' @param alpha Weight parameter: a vector of nonnegative entries the same length as `x`.
#' @param log Logical; whether to evaluate the density on the log scale.
#'
#' @return The density or log-density evaluated at the inputs (scalar).
ddirichlet <- function(x, alpha, log = FALSE) {
  ld <- sum(lgamma(alpha)) - lgamma(sum(alpha))
  ld <- ld + sum((alpha-1) * log(x))
  if(!log) 
    ld <- exp(ld)
  
  ld
}


#' Random sampling from the Dirichlet distribution.
#'
#'
#' @param n Number of random draws.
#' @param alpha Weight parameter: a vector of nonnegative entries.
#' @return A matrix of size `n x length(alpha)` of which each row is a random draw.
rdirichlet <- function(n, alpha) {
  K <- length(alpha) # number of categories
  X <- matrix(rgamma(n*K, shape = alpha), K, n)
  
  drop(t(sweep(X, 2, colSums(X), "/")))
}


#' Solve method for variance matrices.
#'
#'
#' Calculates `y = V^{-1} x` where `V` is a variance matrix.
#'
#' @param V Variance matrix.
#' @param x Optional vector or matrix for which to solve system of equations.  If missing calculates inverse matrix.
#' @param ldV Optionally compute log determinant as well.
#' 
#' @return Either a matrix or vector `y = V^{-1} x`, or if `ldV = TRUE`, a list with elements `y` and `ldV = log(det(V))`.
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


#' Checks whether a specific sample of Z data is Proper 
#'
#' 
#' Ensures that all K groups are represented in the sample Z. 
#' 
#' @param Z A length- `N` vector of integers between 1 and `K`
#' @param k An integer value beween 1 and `K` to indicate which Sigma to condition on 
#' 
#' @return bool indicating whether a sample is proper or not
proper_sample<-function(Z, K){
  #initialize result to true 
  result<-TRUE
  
  #if any group is not represented change result to false and break out of loop
  for (i in 1:K){
    if (sum(Z == i) == 0){
      result<- FALSE
      break
    }
  }
  
  result 
}


#' Simulates Testing Data for Sampler 
#'
#' Follows the data generation process assuming that all parameters are known 
#'
#'@param K scalar value indicating total number of groups 
#'@param N scalar value indicating total number of data samples 
#'@param p scalar value indicating dimension size of each data sample 
#'
#'@return List of all simulated data 
sim_data<-function(K, N, p){
  K<-K
  N<-N
  p<-p 
  
  #set virtual counts to p+2 for all groups 
  vk<- rep(p+2, K)
  
  #initialize rho vector to 0 
  rho<-matrix(rep(0, K), nrow=K, ncol=1)
  rho[,1]<-runif(K)
  rho[,1]<-rho[,1]/sum(rho[,1])
  
  #Sample Z from rho 
  Z <- rep(0, N)
  while(TRUE){
    for (i in 1:N){
      Z[i]<-rcategorical(rho)
    }
    if (proper_sample(Z, K)) break
  } #ensures that sample is proper 
  
  #Simulate Sigma Matrix for all K groups, and ensuring it is positive semi-definite. 
  Sigma <- array(runif(p*p*K)*2-1, dim=c(p,p,K))
  for ( i in 1:K){
    Sigma[,,i]<- t(Sigma[,,i]) %*% Sigma[,,i]
  }
  
  #Simulate V Matrix for all N data points, and ensuring it is positive semi-definite. 
  V <- array(runif(p*p*N)*2-1, dim=c(p,p,N))
  for ( i in 1:N){
    V[,,i]<- t(V[,,i]) %*% V[,,i]
  }
  #Simulate mu values for all K groups
  mu<- array(runif(K*p), dim =c(K, p))
  
  
  #Initialize theta values to all 0
  theta<- matrix( runif(0, p*N), nrow =p, ncol=N)
  
  #Simulate theta values for all N data points by conisdering 
  #their group membership and sampling using the mean and Sigma of that group
  for(i in 1:N){
    z_i <- Z[i]
    theta[,i] <- rmNorm(n=1, mu[z_i,], Sigma[,,z_i])
  }
  
  #Initialize y values (data) to all 0
  y<- matrix( rep(0, p*N), nrow =N, ncol=p)
  
  #Simulate values for all N data points by sampling using their 
  #their corresponding theta and V values. 
  for(i in 1:N){
    y[i,]<- rmNorm(n=1, theta[,i], V[,,i])
  }
  
  #Simulate Omega by using the Variance of the simulated data, y 
  #Omega is the same for all groups. 
  Omega<- array(rep(0,p*p*K), dim=c(p,p,K))
  vary<- var(y)
  
  for(i in 1:K)
    Omega[,,i]<-vary
  
  #Prepare data in a list format 
  data = list(rho, y , V, Z, theta, mu, Sigma, vk, Omega)
  names(data)<-c('rho', 'y', 'V', 'Z', 'theta', 'mu', 'Sigma', 'vk', 'Omega')
  
  data
}


#' @param Z A length- `N` vectir of integers between 1 and `K`
#' @param rho A length-`K` probability vector of group membership.
#' @param theta A `p x N` vector of observation means.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param rand_idx A scalar value bewteen `1` and `N` to indicate which z to condtion on. 
#' @param log bool value indicating whether to give log value of distribution 
zFull.post<-function(Z, rho, theta, mu, Sigma, rand_idx, log=TRUE){
  K<-length(rho)
  
  ll<- rep(0, K)
  ll<- ll + log(rho)
  
  theta_i<-theta[,rand_idx]
  
  for(kk in 1:K){
    ll[kk]<- ll[kk] + dmNorm(theta_i, mu= mu[kk,], Sigma= Sigma[,,kk], log=TRUE)
  }
  
  if (!log)
    ll<-exp(ll)
  
  ll
}


#' @param y
#' @param iter
#' @param init_vals
#' @param burnin_period
#' @param mu.fixed
#' @param Sigma.fixed
#' @param theta.fixed
#' @param rho.fixed
#' @param Z.fixed
#' @param z_out
gibbsSamples <- function(y, iter, init_vals, burnin_period=0, mu.fixed = FALSE, 
                         Sigma.fixed =FALSE, theta.fixed= FALSE, rho.fixed = FALSE, Z.fixed = FALSE, z_out = FALSE ){
  N<-nrow(y)
  p<-ncol(y)

  if(missing(burnin_period)) burnin_period <- min(floor(iter/10), 1e3)
  print(burnin_period)
  mu <- init_vals$mu  
  K<- nrow(mu)
  Sigma <-init_vals$Sigma 
  rho <-init_vals$rho  
  V <-init_vals$V 
  theta <-init_vals$theta 
  Z <-init_vals$Z
  vk<-init_vals$vk 
  Omega<-init_vals$Omega

  mu_samples<- array(NA, dim= c(iter, K, p))
  Sigma_samples<- array(NA, dim = c(iter, p,p, K))
  rho_samples<- array(NA, dim= c(iter, K))
  theta_samples<-array(NA, dim = c(iter, p, N))
  
  if (z_out){
    Z_samples<- matrix (NA, dim = c(iter, N) )
  }
  
  lambda_est <- matrix(rep(0, N*K), nrow=N, ncol=K)
  stopifnot(proper_sample(Z, K))
  
  for(i in 1:(iter+burnin_period)){
    
    if ( mod(i, 50) == 0) {
      cat("Iteration #", i, "\n")
    }
    
    if (!theta.fixed){
      #print("Entering Theta\n")
      theta_Sigma<-array(0, dim=c(p,p,N))
      theta_mu<-array(0, dim=c(N,p))
      for(j in 1:N){
         z_i <- Z[j]
         theta_mu[j,] <- mu[z_i, ]
         theta_Sigma[,,j] <- Sigma[,,z_i]
      }
      theta <-rRxNorm(N, y, V, theta_mu, theta_Sigma)
      theta<-t(theta)
      #print("Exiting Theta\n")
    }
    ###############################################################################
    if (!mu.fixed){
      #print("Entering mu \n")
      Nk<-rep(0, K)
      theta_avg_k<- matrix(rep(0, K*p), nrow=K, ncol=p)
      Sigma_cpy <-Sigma 
      
      for(j in 1:K){
        tidx<- (Z== j)
        Nk[j]<- sum(tidx)
        if (Nk[j] == 1){
          theta_avg_k[j,]<- theta[,tidx]  
        }
        else{
          theta_avg_k[j,]<-rowMeans(theta[,tidx])
        }
        
        Sigma_cpy[,,j]<- Sigma_cpy[,,j]/Nk[j]
      }
      
      mu<-rmNorm(K, theta_avg_k, Sigma_cpy)
      #print("Exiting mu \n")
    }
    ###############################################################################
    if (!Sigma.fixed){
      #print("Entering Sigma \n")
      outer_prod<-array(rep(0, K*p*p), dim=c(p,p,K))
      for(j in 1:K){
        tidx<- (j == Z)
        cur_theta<- theta
        cur_theta_less_mu<- cur_theta[,tidx] - mu[j,]
        outer_prod[,,j]<- cur_theta_less_mu %*% t(cur_theta_less_mu) + Omega[,,j]
      }
      Sigma<-riwish(K, outer_prod, Nk + vk)
      #print("Exiting Sigma \n")
    }
    ################################################################################
    if (!rho.fixed){
      #print("Entering rho \n")
      alpha<- (Nk+1) 
      rho<-rdirichlet(1, alpha)
      #print("Exiting rho \n")
    }
    ###############################################################################
    if(!Z.fixed){
      #print("Entering Z \n")
      kappa<-matrix(rep(0, K*N), nrow=N, ncol=K)
  
      for (n in 1:N){
        kappa[n,]<- zFull.post(Z, rho, theta, mu, Sigma, n)
      }
      
      # if (any(is.na(kappa))){
      #   stop('kappa'); 
      # }
      
      lambda<- exp(kappa)
      
      # if (any(is.na(lambda))){
      #   print(lambda)
      #   stop('lambda'); 
      # }
    
      for(n in 1:N){
        total_exp <- sum(lambda[n,])
        lambda[n,] = lambda[n,]/ total_exp
      }
      
      while(TRUE){
        Z<- rcategorical(t(lambda))
        if (proper_sample(Z, K)) break 
        #print("RESAMPLING Z!")
      }
      
      if (i > burnin_period){
        new_idx <- i-burnin_period
        lambda_est<- (lambda_est*(new_idx -1) + lambda)/new_idx
      }
      #print("Exiting Z \n")
    }
    ###############################################################################
    
    if (i > burnin_period){
      new_idx <- i-burnin_period
      mu_samples[new_idx,,]<-mu
      Sigma_samples[new_idx,,,]<- Sigma
      rho_samples[new_idx,]<- rho
      theta_samples[new_idx,,]<- theta
      
      if (mod(new_idx, 500) == 0){
        print("Mu:")
        print(mu_samples[new_idx,,])
        print("Sigma:")
        print(Sigma_samples[new_idx,,,])
        print("Rho:")
        print(rho_samples[new_idx,])
        print("Theta:")
        print(theta_samples[new_idx,,])
      }
      
      if (z_out)
        Z_samples[new_idx,] <- Z
    }
  }
  
  final_output<-list(mu= mu_samples,  Sigma= Sigma_samples, theta= theta_samples, rho = rho_samples, lambda = lambda_est)
  
  if (z_out)
    final_output<-c(final_output, list(Z = Z_samples))
  
  final_output
}

