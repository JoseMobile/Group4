require(mniw)
require(numbers)
#' Complete data loglikelihood for the normal-normal-mixture (NNM) model.
#'
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param rho A length-`K` probability vector of group membership.
#' @param y An `N x p` matrix of observations, each of which is a column.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#' @param theta A `p x N` matrix of observation means.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param log bool value indicating whether to give log value of distribution 
nnm_loglik <- function(mu, Sigma, rho, y, V, theta, z, log=TRUE) {
  K <- nrow(mu) # number of groups
  N <- nrow(y)
  ll <- 0
  for(kk in 1:K) {
    idk <- z == kk # T/F membership in group k
    Nk <- sum(idk) # number of observations in group k
    #trueIdx<- which(idk == TRUE)
    if(Nk > 0) {
        # group membership contribution
        ll <- ll + Nk * log(rho[kk])
        # mean contribution
        ll <- ll + sum(dmNorm(t(theta[,idk]), mu = matrix(rep(mu[kk,], Nk), nrow=Nk, ncol=p, byrow=TRUE), Sigma = array(rep(Sigma[,,kk], Nk), dim=c(p,p,Nk)),  log = TRUE ))
    }
  }
  # observation contribution
  ll<- ll + sum(dmNorm(y, mu = t(theta), Sigma = V, log = TRUE))
  
  if (!log) {
    ll<- exp(ll)
  }
  
  ll
}

#' Complete data loglikelihood for the normal-normal-mixture (NNM) model.
#'
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param rho A length-`K` probability vector of group membership.
#' @param y An `N x p` matrix of observations, each of which is a column.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#' @param theta A `p x N` matrix of observation means.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param vk A length- `K` vector of integers acting as virtual counts from prior 
#' @param log bool value indicating whether to give log value of distribution 
nnm_post <- function(mu, Sigma, rho, y, V, theta, z, vk, log=TRUE) {
  K<- nrow(mu)
  ll<- nnm_loglik(mu, Sigma, rho, y, V, theta, z)
  for(i in 1:K){
    ll<- ll - 0.5* vk[i]*log(det(Sigma[,,i]))
  }
  if (!log) {ll <- exp(ll)}
  
  ll 
}

#normalized conditional log-pdf
#' @param theta A `p x N` matrix of observation means.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param k An integer value beween 1 and `K` to ondicate which mu to condition on 
#' @param log bool value indicating whether to give log value of distribution 
mu_k.post<-function(theta, mu, Sigma, z, k, log=TRUE){
  mu_k <- mu[k,]
  Sigma_k <- Sigma[,,k]
  tidx <-  tidx<-(k == z)
  Nk<- sum(tidx)
  theta_avg_k<- theta[,tidx]  
  
  if (Nk > 1){
    theta_avg_k<-rowMeans(theta[,tidx])
  }
  
  theta_avg_k<-as.vector(theta_avg_k)
  #ll<-sum(dmNorm(t(theta[,tidx]), log = TRUE,
   #              mu = matrix(rep(mu_k, Nk), nrow=Nk, ncol=p, byrow=TRUE), Sigma = array(rep(Sigma[,,k], Nk), dim=c(p,p,Nk)) ))
  
  ll<-dmNorm(mu_k, theta_avg_k, Sigma_k/Nk, log=TRUE)
  
  if(!log){
    ll<-exp(ll)
  }
  
  ll
}

#' @param theta A `p x N` matrix of observation means.
#' @param marginal_idx 
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param k An integer value beween 1 and `K` to ondicate which mu to condition on 
mu_k_marginals.param<-function(theta, marginal_idx, Sigma, z, k){
  
  Sigma_k <- Sigma[,,k]
  tidx <-  tidx<-(k == z)
  Nk<- sum(tidx)
  theta_avg_k<- theta[,tidx]  
  
  if (Nk > 1){
    theta_avg_k<-rowMeans(theta[,tidx])
  }
  
  theta_avg_k<-as.vector(theta_avg_k)
  Sigma_k <- Sigma_k/Nk
  
  params<- list(mean= theta_avg_k[marginal_idx], sd=sqrt(Sigma_k[marginal_idx, marginal_idx]))
  params
}

#' @param y An `N x p` vector of observations, each of which is a column.
#' @param theta A `p x N` vector of observation means.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param Z A length- `N` vector of integers between 1 and `K`
#' @param rand_idx scalar value netween `1` and `N` to indicate which theta to condition on 
#' @param log bool value indicating whether to give log value of distribution 
theta_i.post<-function(y, theta, V, mu, Sigma, Z, rand_idx, log=TRUE){
  y_i <- y[rand_idx,]
  theta_i<- theta[,rand_idx]
  V_i<- V[,, rand_idx]
  z_i <- Z[rand_idx]
  Sigma_zi <- Sigma[,,z_i]
  mu_zi <- mu[z_i,]
  
  #ll<- dmNorm(theta_i,  mu = mu_zi, Sigma = Sigma_zi, log = TRUE ) 
  #ll<- ll + dmNorm(y_i, mu = theta_i, Sigma = V_i, log=TRUE )
  
  G_i <- Sigma_zi %*% solveV(V_i + Sigma_zi)
  new_mu <- G_i  %*% (y_i - mu_zi) + mu_zi
  new_V <- G_i  %*% V_i 
  
  new_mu<-as.vector(new_mu)
  ll<- dmNorm(theta_i, new_mu, new_V, log=TRUE)
  
  if(!log){
    ll<-exp(ll)
  }
  
  ll
}

#' @param theta A `p x N` matrix of observation means.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param Z A length- `N` vector of integers between 1 and `K`
#' @param k An integer value beween 1 and `K` to indicate which Sigma to condition on 
#' @param vk A length- `K` vector of integers acting as virtual counts from prior 
#' @param log bool value indicating whether to give log value of distribution 
sigma_k.post<-function(theta, mu, Sigma, Z, k, vk, log=TRUE){
  Sigma_k <- Sigma[,, k]
  Sigma_inv<- solveV(Sigma_k)
  mu_k <- mu[k, ]
  vkk<- vk[k]
  tidx<- (Z == k)
  Nk<- sum(tidx)
  theta_less_mu <- theta[,tidx] - mu_k
  ll <- diwish(Sigma_k, (theta_less_mu %*% t(theta_less_mu)) ,  Nk + vkk -p -1, log = TRUE)
  
  #ll<- ll - 0.5* sum(diag(Sigma_inv %*% (theta_less_mu %*% t(theta_less_mu))))
  #ll <- ll - 0.5*(Nk + vkk)*log(det(Sigma_k))
  
  if (!log){
    ll<-exp(ll)
  } 
  
  ll
}

#' @param Z A length- `N` vectir of integers between 1 and `K`
#' @param rho A length-`K` probability vector of group membership.
#' @param log bool value indicating whether to give log value of distribution 
rho.post<-function(Z, rho, log=TRUE){
  K<- length(rho)
  Nk<- rep(0, K)
  log_rho <- log(rho)
  for(j in 1:K){
    tidx<- (j == Z)
    Nk[j] <- sum(tidx)
  }
  #ll<-sum(Nk*log_rho)
  ll<-ddirichlet(rho, Nk + 1, log=TRUE)
  
  if (!log)
    ll<-exp(ll)
  
  ll
}

#' @param Z A length- `N` vectir of integers between 1 and `K`
#' @param rho A length-`K` probability vector of group membership.
#' @param theta A `p x N` vector of observation means.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param rand_idx A scalar value bewteen `1` and `N` to indicate which z to condtion on. 
#' @param log bool value indicating whether to give log value of distribution 
z.post<-function(Z, rho, theta, mu, Sigma, rand_idx, log=TRUE){
  ll<-0 
  z_i <- Z[rand_idx]
  theta_i <- theta[, rand_idx]
  ll<- ll + log(rho[z_i])
  
  Sigma_k <- Sigma[,, z_i]
  mu_k <- mu[z_i,]
  
  ll<- ll + dmNorm(theta_i, mu= mu_k, Sigma= Sigma_k, log=TRUE)
  
  if (!log)
    ll<-exp(ll)
  
  ll
}


#' Sample from a categorical distribution.
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
#' @param x Observation vector of nonnegative entries which sum to one.
#' @param alpha Weight parameter: a vector of nonnegative entries the same length as `x`.
#' @param log Logical; whether to evaluate the density on the log scale.
#'
#' @return The density or log-density evaluated at the inputs (scalar).
ddirichlet <- function(x, alpha, log = FALSE) {
  ld <- sum(lgamma(alpha)) - lgamma(sum(alpha))
  ld <- ld + sum((alpha-1) * log(x))
  if(!log) ld <- exp(ld)
  ld
}

#' Random sampling from the Dirichlet distribution.
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
#' Calculates `y = V^{-1} x` where `V` is a variance matrix.
#'
#' @param V Variance matrix.
#' @param x Optional vector or matrix for which to solve system of equations.  If missing calculates inverse matrix.
#' @param ldV Optionally compute log determinant as well.
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

#' @param Z A length- `N` vector of integers between 1 and `K`
#' @param k An integer value beween 1 and `K` to indicate which Sigma to condition on 
proper_sample<-function(Z, K){
  result<-TRUE
  for (i in 1:K){
    if (sum(Z == i) == 0){
      result<- FALSE
    }
  }
  result 
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
  mu <- init_vals$mu  
  K<- nrow(mu)
  Sigma <-init_vals$Sigma 
  rho <-init_vals$rho  
  V <-init_vals$V 
  theta -init_vals$theta 
  Z <-init_vals$Z
  vk<init_vals$vk 
  
  mu_samples<- array(NA, dim= c(iter, K, p))
  Sigma_samples<- array(NA, dim = c(iter, p,p, K))
  rho_samples<- array(NA, dim= c(iter, N))
  theta_samples<-array(NA, dim = c(iter, p, N))
  
  if (z_out){
    Z_samples<- matrix (NA, dim = c(iter, N) )
  }
  
  lambda_est <- matrix(rep(0, N*K), nrow=N, ncol=K)
  stopifnot(proper_sample(Z, K))
  
  for(i in 1:(iter+burnin_period)){
    
    if ( mod(i, 500) == 0) {
      cat("Iteration #", i, "\n")
    }
    
    if (!theta.fixed){
      theta <-rRxNorm(N, y, V, mu, Sigma)
    }
    ###############################################################################
    if (!mu.fixed){
      Nk<-rep(0, K)
      theta_avg_k<- matrix(rep(0, K*p), nrow=K, ncol=p)
      Sigma_cpy <-Sigma 
      
      for(j in 1:K){
        tidx<- (j == Z)
        Nk[j]<- sum(tidx)
        if (Nk[j] == 1)
          theta_avg_k[j,]<- theta[,tidx]  
        else
          theta_avg_k[j,]<-rowMeans(theta[,tidx])
        
        Sigma_cpy[,,j]<- Sigma_cpy[,,j]/Nk[j]
      }
      
      mu<-rmNorm(K, theta_avg_k, Sigma_cpy)
    }
    ###############################################################################
    if (!Sigma.fixed){
      outer_prod<-array(rep(0, K*p*p), dim=c(p,p,K))
      for(j in 1:K){
        tidx<- (j == Z)
        cur_theta<-theta
        cur_theta_less_mu<- cur_theta[,tidx] - mu[j,]
        outer_prod[,,i]<- cur_theta_less_mu %*% t(cur_theta_less_mu)
      }
      Sigma<-riwish(K, outer_prod, Nk + vk -p-1)
    }
    ################################################################################
    if (!rho.fixed){
      alpha<-Nk+1 
      rho<-rdirichlet(1, alpha)
    }
    ###############################################################################
    if(!Z.fixed){
      kappa<-matrix(rep(0, K*N), nrow=N, ncol=K)
  
      for(j in 1:K){
        for (n in 1:N){
          kappa[n,j]<- z.post(Z, rho, theta, mu, Sigma, n)
        }
      }
      
      lambda<- exp(kappa)
      
      for(n in 1:N){
        total_exp <- sum(lambda[n,])
        lambda[n,] = lambda[n,]/ total_exp
      }
      
      while(TRUE){
        Z<- rcategorical(t(lambda))
        if (proper_sample(Z, K)) break 
      }
      
      if (i > burnin_period)
        lambda_est<- (lambda_est*(new_idx -1) + lambda)/new_idx
    }
    ###############################################################################
    
    if (i > burnin_period){
      new_idx <- i-burnin_period
      mu_samples[new_idx,,]<-mu 
      Sigma_samples[new_idx,,,]<- Sigma
      rho_samples[new_idx,]<- rho
      theta_samples[new_idx,,]<- theta
      
      if (z_out)
        Z_samples[new_idx,] <- Z
    }
  }
  
  final_output<-list(mu= mu_samples,  Sigma= Sigma_samples, theta= theta_samples, rho = rho_samples, lambda = lambda_est)
  
  if (z_out)
    final_output<-c(final_output, list(Z = Z_samples))
  
  final_output
}

