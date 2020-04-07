require(mniw)

#' Complete data loglikelihood for the normal-normal-mixture (NNM) model.
#'
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param rho A length-`K` probability vector of group membership.
#' @param y An `N x p` matrix of observations, each of which is a column.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#' @param theta A `p x N` matrix of observation means.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
nnm_loglik <- function(mu, Sigma, rho, y, V, theta, z) {
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
  ll + sum(dmNorm(y, mu = t(theta), Sigma = V, log = TRUE))
}

#normalized conditional log-pdf
#' @param theta A `p x N` matrix of observation means.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param k An integer value beween 1 and `K` to ondicate which mu to condition on 
mu_k.loglik<-function(theta, mu, Sigma, z, k){
  mu_k <- mu[k,]
  tidx <-  tidx<-(k == z)
  Nk<- sum(tidx)
  sum(dmNorm(t(theta[,tidx]), log = TRUE,
                 mu = matrix(rep(mu_k, Nk), nrow=Nk, ncol=p, byrow=TRUE), Sigma = array(rep(Sigma[,,k], Nk), dim=c(p,p,Nk)) ))
}

#' @param y An `N x p` vector of observations, each of which is a column.
#' @param theta A `p x N` vector of observation means.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param Z A length- `N` vector of integers between 1 and `K`
#' @param rand_idx scalar value netween `1` and `N` to indicate which theta to condition on 
theta_i.loglik<-function(y, theta, V, mu, Sigma, Z, rand_idx){
  y_i <- y[rand_idx,]
  theta_i<- theta[,rand_idx]
  V_i<- V[,, rand_idx]
  z_i <- Z[rand_idx]
  Sigma_zi <- Sigma[,,z_i]
  mu_zi <- mu[z_i,]
  
  ll<- 0 
  ll<- ll + dmNorm(theta_i,  mu = mu_zi, Sigma = Sigma_zi, log = TRUE )
  ll<- ll  + dmNorm(y_i, mu = theta_i, Sigma = V_i, log=TRUE )
  ll
}

#' @param theta A `p x N` matrix of observation means.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param Z A length- `N` vectir of integers between 1 and `K`
#' @param k An integer value beween 1 and `K` to indicate which Sigma to condition on 
sigma_k.loglik<-function(theta, mu, Sigma, Z, k){
  ll<- 0 
  Sigma_k <- Sigma[,, k]
  Sigma_inv<- solveV(Sigma_k)
  mu_k <- mu[k, ]
  tidx<- (Z == k)
  Nk<- sum(tidx)
  theta_less_mu <- theta[,tidx] - mu_k
  ll<- ll - 0.5* sum(diag(Sigma_inv %*% (theta_less_mu %*% t(theta_less_mu))))
  ll <- ll - 0.5*Nk*log(det(Sigma_k))
  ll
}

#' @param Z A length- `N` vectir of integers between 1 and `K`
#' @param rho A length-`K` probability vector of group membership.
rho.loglik<-function(Z, rho){
  K<- length(rho)
  Nk<- rep(0, K)
  log_rho <- log(rho)
  for(j in 1:K){
    tidx<- (j == Z)
    Nk[j] <- sum(tidx)
  }
  sum(Nk*log_rho)
}

#' @param Z A length- `N` vectir of integers between 1 and `K`
#' @param rho A length-`K` probability vector of group membership.
#' @param theta A `p x N` vector of observation means.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param rand_idx A scalar value bewteen `1` and `N` to indicate which z to condtion on. 
z.loglik<-function(Z, rho, theta, mu, Sigma, rand_idx){
  ll<-0 
  z_i <- Z[rand_idx]
  theta_i <- theta[, rand_idx]
  ll<- ll + log(rho[z_i])
  
  Sigma_k <- Sigma[,, z_i]
  mu_k <- mu[z_i,]
  
  ll<- ll + dmNorm(theta_i, mu= mu_k, Sigma= Sigma_k, log=TRUE)
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

#unfinished Sampler 
gibbsSamples <- function(y, K, iter, init_vals ){
  N<-nrow(y)
  p<-ncol(y)
  
  mu_Iter<-array(rep(0, iter*K*p), dim=c(iter, K, p))
  mu_Iter[1,,]<- init_vals$mu 
  
  Sigma_Iter<-array(rep(0, iter*K*p*p), dim=c(iter, p, p, K))
  Sigma_Iter[1,,,]<-init_vals$Sigma 
  
  rho_Iter<-array(rep(0, iter*K), dim=c(iter, K))
  rho_Iter[1,]<-init_vals$rho 
  
  V_Iter<-array(rep(0, iter*N*p*p), dim=c(iter, p, p, N))
  V_Iter[1,,,]<-init_vals$V 
  
  theta_Iter<-array(rep(0, p*N*iter), dim=c(iter,p,N))
  theta_Iter[1,,]<-init_vals$theta 
  
  z_Iter<-array(rep(0, N*iter), dim=c(iter,N))
  z_Iter[1,]<-init_vals$z
  
  for( i in 2:iter){
    theta_Iter[i,,]<-rRxNorm(1, y, V_Iter[i-1,,,], mu_Iter[i-1,,], Sigma[i-1,,,])
    ###############################################################################
    Nk<-rep(0, 10)
    theta_avg_k<- matrix(rep(0, K*p), nrow=K, ncol=p)
    for(j in 1:K){
      tidx<- (j == z_Iter[i-1,])
      Nk[j]<- sum(tidx)
      theta_avg_k[j,]<-rowMeans(theta_Iter[i,,tidx])
    }
    mu_Iter[i,,]<-rmNorm(1, theta_avg_k, Sigma[i-1,,,]/Nk)
    ###############################################################################
    Psi<-array(rep(0, K*p*p), dim=c(p,p,K))
    for(j in 1:K){
      tidx<- (j == z_Iter[i-1,])
      cur_theta<-theta_Iter[i,,]
      cur_theta_k<- cur_theta[,tidx] - mu_Iter[i,j,]
      Psi[,,i]<- cur_theta_k %*% t(cur_theta_k)
    }
    Sigma_Iter[i,,,]<-riwish(1, Psi, Nk-p-1)
    ################################################################################
    alpha<-Nk+1 
    rho_Iter[i,]<-rdirichlet(1, alpha)
    ###############################################################################
    kappa<-matrix(rep(0, K*N), nrow=N, ncol=K)
    Sigma_inv<-array(rep(0, K*p*p), dim=c(p, p, K))
    
    for(j in 1:K){
      Sigma_inv[,,j]<- solveV(Sigma_Iter[i,,,j])
      log_rho_k<-log(rho_Iter[i,j])
      log_det_k<-log(det(Psi[,,j]))

      for (n in 1:N){
        theta_sub_mu <- theta_Iter[i,,n] - mu_Iter[i,j,] 
        kappa[n,j]<- log_rho_k - 0.5 *( t(theta_sub_mu) %*% Sigma_inv[,,j] %*% theta_sub_mu + log_det_k)
      }
    }
    
    lambda<- exp(kappa)
    assert(sum(is.na(lambda)) == 0)
    
    for(n in 1:N){
      total_exp <- sum(lambda[n,])
      lambda[n,] = lambda[n,]/ total_exp
    }
    Z_iter[i,] <- rcategorical(t(lambda))
  }
}

