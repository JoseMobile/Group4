require(mniw)

# =============================================================================
# FUNCTIONS FOR LIKELIHOODS AND POSTERIORS
# =============================================================================

#' Density of the Dirichlet distribution.
#'
#'
#' @param x Observation vector of nonnegative entries which sum to one.
#' @param alpha Weight parameter: a vector of nonnegative entries the same length as `x`.
#' @param log Logical; whether to evaluate the density on the log scale.
#'
#' @return The density or log-density evaluated at the inputs (scalar).
#' @author Martin Lysy
ddirichlet <- function(x, alpha, log = FALSE) {
  ld <- sum(lgamma(alpha)) - lgamma(sum(alpha))
  ld <- ld + sum((alpha-1) * log(x))
  if(!log)
    ld <- exp(ld)

  ld
}

#' Complete data loglikelihood for the normal-normal-mixture (NNM) model.
#'
#'
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param rho A length-`K` probability vector of group membership.
#' @param y An `N x p` matrix of observations, each of which is a column.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#' @param theta A `N x p` matrix of observation means.
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
      ll <- ll + sum(dmNorm(theta[idk,], mu = matrix(rep(mu[kk,], Nk), nrow=Nk, ncol=p, byrow=TRUE), Sigma = array(rep(Sigma[,,kk], Nk), dim=c(p,p,Nk)),  log = TRUE ))
    }
  }
  # observation contribution
  ll<- ll + sum(dmNorm(y, mu =theta, Sigma = V, log = TRUE))

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
#' @param theta A `N x p` matrix of observation means.
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
#' @param theta A `N x p` matrix of observation means.
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
  theta_avg_k<- theta[tidx,]

  #Nk > 1, avg group members values across each dimension
  if (Nk > 1){
    theta_avg_k<-colMeans(theta[tidx,])
  }

  #ensure vector data structure is used
  theta_avg_k<-as.vector(theta_avg_k)

  #get log pdf of a multivariate normal distribution
  ll<-dmNorm(mu_k, theta_avg_k, Sigma_k/Nk, log=TRUE)

  if(!log) # check if non log value desired
    ll<-exp(ll)

  ll
}

#' Normalized Conditional Log-PDF of Theta_i
#'
#'
#' @param y An `N x p` vector of observations, each of which is a column.
#' @param theta A `N x p` vector of observation means.
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
  theta_i<- theta[rand_idx,]
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
#' @param theta A `N x p` matrix of observation means.
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

  theta<-t(theta)
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
#' @param theta A `N x p` vector of observation means.
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
  theta_i <- theta[rand_idx,]

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
# =============================================================================

# =============================================================================
# SIMULATED DATA FUNCTION sim_data
# =============================================================================

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
  theta<- matrix( runif(0, p*N), nrow =N, ncol=p)

  #Simulate theta values for all N data points by conisdering
  #their group membership and sampling using the mean and Sigma of that group
  for(i in 1:N){
    z_i <- Z[i]
    theta[i,] <- rmNorm(n=1, mu[z_i,], Sigma[,,z_i])
  }

  #Initialize y values (data) to all 0
  y<- matrix( rep(0, p*N), nrow =N, ncol=p)

  #Simulate values for all N data points by sampling using their
  #their corresponding theta and V values.
  for(i in 1:N){
    y[i,]<- rmNorm(n=1, theta[i,], V[,,i])
  }

  #Simulate Omega by using the Variance of the simulated data, y
  #Omega is the same for all groups.
  Omega<- array(rep(0,p*p*K), dim=c(p,p,K))
  vary<- var(y)*(N-1)/(N-p)

  for(i in 1:K)
    Omega[,,i]<-vary

  #Prepare data in a list format
  data = list(rho, y , V, Z, theta, mu, Sigma, vk, Omega)
  names(data)<-c('rho', 'y', 'V', 'Z', 'theta', 'mu', 'Sigma', 'vk', 'Omega')

  data
}
# =============================================================================

###################################  Testing Mu ################################################
K<-3
N<-100
p<-3
num_sim =10
data<- sim_data(K, N, p)
rho <- data$rho
y <- data$y
V <-data$V
Z <- data$Z
theta <- data$theta
mu <-data$mu
Sigma <-data$Sigma
vk<-data$vk
Omega<-data$Omega

rand_idx <- ceiling(K*runif(1))

ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  mu[rand_idx,]<- runif(p)
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk, Omega)
  lld[i] <- mu_k.post(theta, mu, Sigma, Z, rand_idx)
}

difference <- lld-ldu

# iterate through all values of difference and make sure
# they are equal to the first value
for (diff_i in difference) {
  test_that("mu posterior", {
    expect_equal(diff_i, difference[1])
  })
}

###################################  Testing Theta ################################################

data<- sim_data(K, N, p)
rho <- data$rho
y <- data$y
V <-data$V
Z <- data$Z
theta <- data$theta
mu <-data$mu
Sigma <-data$Sigma
vk<-data$vk
Omega<-data$Omega

rand_idx <- ceiling(N*runif(1))

ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  z_i <- Z[rand_idx]
  theta[rand_idx,]<- rmNorm(1, mu[z_i,], Sigma[,,z_i])
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk, Omega)
  lld[i] <- theta_i.post(y, theta, V, mu, Sigma, Z, rand_idx)
}

difference <- lld-ldu

# iterate through all values of difference and make sure
# they are equal to the first value
for (diff_i in difference) {
  test_that("theta posterior", {
    expect_equal(diff_i, difference[1])
  })
}

###################################  Testing Sigma ################################################

data<- sim_data(K, N, p)
rho <- data$rho
y <- data$y
V <-data$V
Z <- data$Z
theta <- data$theta
mu <-data$mu
Sigma <-data$Sigma
vk<-data$vk
Omega<-data$Omega

rand_idx <- ceiling(K*runif(1))

ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  new_Sigma <-matrix(runif(p*p)*2-1, nrow= p, ncol =p )
  new_Sigma <- t(new_Sigma) %*% new_Sigma
  Sigma[,,rand_idx] <- new_Sigma
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk, Omega)
  lld[i] <- sigma_k.post(theta, mu, Sigma, Z, rand_idx, vk, Omega)
}

difference <- lld-ldu

# iterate through all values of difference and make sure
# they are equal to the first value
for (diff_i in difference) {
  test_that("Sigma posterior", {
    expect_equal(diff_i, difference[1])
  })
}

###################################  Testing rho ################################################

data<- sim_data(K, N, p)
rho <- data$rho
y <- data$y
V <-data$V
Z <- data$Z
theta <- data$theta
mu <-data$mu
Sigma <-data$Sigma
vk<-data$vk
Omega<-data$Omega


ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  rho[,1]<-runif(K)
  rho[,1]<-rho[,1]/sum(rho[,1])
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk, Omega)
  lld[i] <- rho.post(Z, rho[,1])
}

difference <- lld-ldu

# iterate through all values of difference and make sure
# they are equal to the first value
for (diff_i in difference) {
  test_that("rho posterior", {
    expect_equal(diff_i, difference[1])
  })
}

###################################  Testing Z ################################################

data<- sim_data(K, N, p)
rho <- data$rho
y <- data$y
V <-data$V
Z <- data$Z
theta <- data$theta
mu <-data$mu
Sigma <-data$Sigma
vk<-data$vk
Omega<-data$Omega


rand_idx <- ceiling(N*runif(1))

ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  Z[rand_idx]<- rcategorical(rho)
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk, Omega)
  lld[i] <- z.post(Z, rho[,1], theta, mu, Sigma, rand_idx)
}

difference <- lld-ldu

# iterate through all values of difference and make sure
# they are equal to the first value
for (diff_i in difference) {
  test_that("z posterior", {
    expect_equal(diff_i, difference[1])
  })
}

