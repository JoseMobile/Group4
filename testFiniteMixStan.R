
require(rstan)
require(mniw)
source('nnm-functions.R')
source("gibbs-wrapper.R")

hbe_fit <- readRDS("hbe-fbm_fit_v1.rds")
N <- length(hbe_fit)
p <- length(hbe_fit[[1]]$coef)

# Extract data 
wts <- vector(length=N)
theta_hat <- matrix(nrow=N, ncol=p)
obs_info <- array(dim=c(p, p, N))
for (i in 1:N){
  theta_hat[i, ] <- hbe_fit[[i]]$coef
  wts[i] <- hbe_fit[[i]]$wt
  obs_info[, , i] <- hbe_fit[[i]]$vcov
}

# Transform the wts into numeric labels
lab <- unique(wts) 
K <- length(lab)
cluster <- numeric(length=N)
for (i in 1:N){
  if (wts[i] == lab[1]){
    cluster[i] <- 1
  }
  else if (wts[i] == lab[2]){
    cluster[i] <- 2
  }
  else if (wts[i] == lab[3]){
    cluster[i] <- 3
  }
  else if (wts[i] == lab[4]){
    cluster[i] <- 4
  }
  else if (wts[i] == lab[5]){
    cluster[i] <- 5
  }
  else{
    cluster[i] <- 6
  }
}

real_N_k <- numeric(length=K)
for (k in 1:K){
  real_N_k[k] <- sum(cluster == k)
}

finMix_mod <- stan_model("finiteMixtures.stan")

p = 6
N = 668 
K = 6 

initZ <- init_z(theta_hat, 6)
count <- as.numeric(table(initZ))
mu_init <- array(0, dim =c(K, p))

for( i in 1:K){
  tidx<- initZ == i
  if (count[i] ==0)
    mu_init[i,]<-colMeans(theta_hat) 
  else if (count[i] == 1)
    mu_init[i,]<- theta_hat[tidx,]
  else {
    mu_init[i,]<- colMeans(theta_hat[tidx,]) 
  }
}

sigma_init <- array(0, dim =c(K,p,p))
for( i in 1:K){
  tidx <- (initZ == i) 
  if (count[i] < p)
    sigma_init[i, , ]<-var(theta_hat) *(N-1)/(N-p)
  else {
    sigma_init[i, , ]<- var(theta_hat[tidx,]) *(count[i] -1)/(count[i]-p)
  }
}

rho_init <- count/N
theta_init <- theta_hat

#Data Reshape to fit Stan Model data structure shape
reshapeV<-array(0, dim=c(N,p,p))
for(i in 1:N){
  reshapeV[i,,]<-obs_info[,,i]
}

reshapeOmega<-array(0, dim=c(K,p,p))
for(i in 1:K){
  reshapeOmega[i,,]<- var(theta_hat) *(N-1)/(N-p)
}
vK<-rep(p+2, K)

#Creating param list 
param_init = list(sigma = sigma_init, mu= mu_init, theta=theta_init, rho=rho_init)
param_list = list(param_init, param_init, param_init)

finMixData <- list(K =K, N = N, p=p, V = reshapeV, Omega=reshapeOmega, vK = vK, Y=theta_hat )


### THIS WILL TAKE 5 + hours to run
#Perform Sampling 
finMixData_fit <- sampling(finMix_mod, data = finMixData, iter = 10000, warmup=2000,
                           verbose = TRUE, chains = 3, init=param_list, control = list(max_treedepth = 10))

saveRDS(finMixData_fit, "fitMixDataStan.rds")

#################################  Checking for Convergence of Stan Model ############################################

#finMixData_fit<- readRDS("fitMixDataStan.rds")
  
print(finMixData_fit)

#Checking Stan Trace 
stan_trace(finMixData_fit, 'rho[3]') 
stan_trace(finMixData_fit, 'sigma[1,2,1]') 
stan_trace(finMixData_fit, 'sigma[2,4,3]') 
stan_trace(finMixData_fit, 'sigma[5,5,5]') 
stan_trace(finMixData_fit, 'mu[2,5]') 
stan_trace(finMixData_fit, 'mu[4,4]') 


#Checking Stan Density
rstan::stan_dens(finMixData_fit, 'rho[3]', separate_chains=TRUE)
rstan::stan_dens(finMixData_fit, 'sigma[1,2,1]', separate_chains=TRUE)
rstan::stan_dens(finMixData_fit, 'sigma[2,4,3]', separate_chains=TRUE)
rstan::stan_dens(finMixData_fit, 'sigma[5,5,5]', separate_chains=TRUE)
rstan::stan_dens(finMixData_fit, 'mu[2,5]', separate_chains=TRUE)
rstan::stan_dens(finMixData_fit, 'mu[4,4]', separate_chains=TRUE)

#Checking ESS Bulk
mu25_sim<-array(extract(finMixData_fit)$mu[,2,5], dim=c(8000,3))
sigma121_sim<-array(extract(finMixData_fit)$sigma[,1,2,1], dim=c(8000,3))
sigma243_sim<-array(extract(finMixData_fit)$sigma[,2,4,3], dim=c(8000,3))
sigma555_sim<-array(extract(finMixData_fit)$sigma[,5,5,5], dim=c(8000,3))
rho3_sim <- array(extract(finMixData_fit)$rho[,3], dim=c(8000,3))
mu44_sim<-array(extract(finMixData_fit)$mu[,4,4], dim=c(8000,3))

ess_bulk(mu25_sim)
ess_bulk(sigma121_sim)
ess_bulk(sigma243_sim)
ess_bulk(sigma555_sim)
ess_bulk(rho3_sim)
ess_bulk(mu44_sim)


#pair plot
pairs(finMixData_fit, pars = c("rho[1]", "rho[2]", "rho[3]", "rho[4]"))

#summary plot
summary(finMixData_fit)$summary


################################# Get cluster groups ###############################################

extract_rho<- extract(finMixData_fit)$rho
extract_theta<-  extract(finMixData_fit)$theta
extract_mu<-  extract(finMixData_fit)$mu
extract_sigma<-  extract(finMixData_fit)$sigma


print(dim(extract_rho))
print(dim(extract_theta))
print(dim(extract_mu))
print(dim(extract_sigma))


numIter<- 24000
Lambda<-matrix(0, nrow =N, ncol=K)


for ( i in 1:numIter){
  cur_rho<- extract_rho[i,]
  cur_mu<- extract_mu[i, ,]
  cur_theta<- extract_theta[i, ,]
  cur_sigma<- extract_sigma[i, , ,]
  K<-length(cur_rho)
  Kappa<-matrix(0, nrow=N, ncol=K)
  inv_sigma <- vector("list", length=K)
  
  for (k in 1:K){
    inv_sigma[[k]] <- solveV(cur_sigma[k,,], ldV=TRUE)
  }
  
  log_det <- sapply(inv_sigma, FUN=function(z){ z$ldV })
  log_rho <- log(cur_rho)
  
  for (j in 1:N){
    for (k in 1:K){
      temp <- as.matrix(cur_theta[j, ] - cur_mu[k, ])
      Kappa[j, k] <- log_rho[k] - (t(temp) %*% inv_sigma[[k]]$y %*% temp + log_det[k]) / 2
    }
  }
    
  cur_Lambda<- exp(Kappa)
  
  if ( sum(is.nan(cur_Lambda)) > 0){
    print("NaN found.")
    stop() 
  }
  
  cur_Lambda <- cur_Lambda / rowSums(cur_Lambda)
  
  Lambda <- Lambda + cur_Lambda
  
  if(mod(i,120) == 0) print(i)
}

Lambda <- Lambda/ numIter

fitted_clusters <- numeric(length=N)
for (i in 1:N){
  fitted_clusters[i] <- which.max(Lambda[i, ])
}

print(table(fitted_clusters))

################################# Sampler Testing  ###############################################

#Begin Testing 
num_sim<-10
fullPost<-rep(0, num_sim)
stanPost<-rep(0, num_sim)



#' Log Sum Exp 
#' More stable way to add a vector of loglilikelihood
#'
#'
#'@param x vector of loglikelihoods to be added 
#'@return scalar sum of loglikelihoods
log_sum_exp <- function(x) {
  xmax <- max(x)
  log(sum(exp(x - xmax))) + xmax
}

#' New Complete Data Posterior of Multivariate Mixture Model (for RStan Comparison)
#' 
#' 
#' @param Y An `N x p` matrix of observations, each of which is a column.
#' @param vK A length- `K` vector of integers acting as virtual counts from prior 
#' @param Omega A `p x p x K` matrix of variances corresponding to the variance of the prior chosen.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#' @param rho A length-`K` probability vector of group membership.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param theta A `N x p` matrix of observation means.
#' @return full log data posterior of model 
nnm2_post<-function(Y, vK, Omega, V, rho, mu, Sigma, theta){
  ll<-0
  ll<- ll + sum(diwish(Sigma, Omega , vK, log = TRUE))
  log_rho<-log(rho)
  for(i in 1:N){
    ll_idv<-log_rho
    theta_i<- theta[i,]
    repTheta_i<-matrix(0, nrow=K, ncol=p)
    for( j in 1:K){
      repTheta_i[j,]<- theta_i
    }
    ll_idv<- ll_idv + dmNorm(repTheta_i, mu, Sigma, log=TRUE)
    ll<- ll + log_sum_exp(ll_idv)
  }
  ll<- ll + sum(dmNorm(Y, theta, V, log=TRUE))
  
  ll
}

#' Get a new set of starting param values
#' 
#' 
#' @param theta_hat  An `N x p` matrix of observations, each of which is a column.
#' @return list of param initialization values using a different starting set of z's(from km++ algorithm)
new_params<-function(theta_hat){
  initZ <- init_z(theta_hat, 6)
  count <- as.numeric(table(initZ))
  mu_init <- array(0, dim =c(K, p))
  
  print(count)
  
  for( i in 1:K){
    tidx<- initZ == i
    if (count[i] ==0)
      mu_init[i,]<-colMeans(theta_hat) 
    else if (count[i] == 1)
      mu_init[i,]<- theta_hat[tidx,]
    else {
      mu_init[i,]<- colMeans(theta_hat[tidx,]) 
    }
  }
  
  sigma_init <- array(0, dim =c(K,p,p))
  for( i in 1:K){
    tidx <- (initZ == i) 
    if (count[i] < p)
      sigma_init[i, , ]<-var(theta_hat) *(N-1)/(N-p)
    else {
      sigma_init[i, , ]<- var(theta_hat[tidx,]) *(count[i] -1)/(count[i]-p)
    }
  }
  
  rho_init <- count/N
  theta_init <- theta_hat
  
  list(theta = theta_init, rho= rho_init, sigma=sigma_init, mu=mu_init)
}


# Perform Simulations 
for(i in 1:num_sim){
  
  #Create New Param Values (We have 4 params)
  new_data<- new_params(theta_hat)
  new_rho <- new_data$rho
  new_theta <- new_data$theta
  new_mu <- new_data$mu
  new_sigma <- new_data$sigma
  
  Omega<-array(0, dim=c(p,p,K))
  for(j in j:K){
    Omega[,,j]<-reshapeOmega[j,,]
  }
  
  rSigma<-array(0, dim=c(p,p,K))
  for(j in j:K){
    rSigma[,,j]<-new_sigma[j,,]
  }
  
  fullPost[i]<- nnm2_post(theta_hat, vK, Omega, obs_info, new_rho, new_mu, rSigma, new_theta)
  ########################################################################
  pars<- list(sigma = new_sigma, mu= new_mu, theta=new_theta, rho=new_rho)
  upars <- unconstrain_pars(finMixData_fit, pars)
  stanPost[i]<-log_prob(object = finMixData_fit, upars = upars, adjust_transform = FALSE)
}

print(fullPost)
print(stanPost)
fullPost-stanPost

