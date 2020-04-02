
require(rstan)
require(mniw)
require(condMVNorm)
source('nnm-functions.R')

###################################  Simulate Data ################################################
#Helper Function to Simulate the Data 

proper_sample<-function(Z, K){
  result<-TRUE
  for (i in 1:K){
    if (sum(Z == i) == 0){
      result<- FALSE
    }
  }
  result 
}

sim_data<-function(K, N, p){
  K<-K
  N<-N
  p<-p 
  Z <- rep(0, N)
  
  vk<- sample( (2*p-1): (4*p- 1), K, replace=TRUE)
  rho<-matrix(rep(0, K), nrow=K, ncol=1)
  rho[,1]<-runif(K)
  rho[,1]<-rho[,1]/sum(rho[,1])
  
  while(TRUE){
    for (i in 1:N){
       Z[i]<-rcategorical(rho)
    }
    if (proper_sample(Z, K)) break
  }
  
  Sigma <- array(runif(p*p*K)*2-1, dim=c(p,p,K))
  for ( i in 1:K){
    Sigma[,,i]<- t(Sigma[,,i]) %*% Sigma[,,i]
  }
  
  V <- array(runif(p*p*N)*2-1, dim=c(p,p,N))
  for ( i in 1:N){
    V[,,i]<- t(V[,,i]) %*% V[,,i]
  }
  
  mu<- array(runif(K*p), dim =c(K, p))
  theta<- matrix( runif(0, p*N), nrow =p, ncol=N)
  
  for(i in 1:N){
    z_i <- Z[i]
    theta[,i] <- rmNorm(n=1, mu[z_i,], Sigma[,,z_i])
  }
  
  y<- matrix( rep(0, p*N), nrow =N, ncol=p)
  for(i in 1:N){
    y[i,]<- rmNorm(n=1, theta[,i], V[,,i])
  }
  
  data = list(rho, y , V, Z, theta, mu, Sigma, vk)
  
  names(data)<-c('rho', 'y', 'V', 'Z', 'theta', 'mu', 'Sigma', 'vk')
  data
}


###################################  Testing Mu ################################################
K<-10
N<-100
p<-3 

data<- sim_data(K, N, p)
rho <- data$rho
y <- data$y
V <-data$V
Z <- data$Z
theta <- data$theta
mu <-data$mu
Sigma <-data$Sigma
vk<-data$vk

rand_idx<- as.integer(K*runif(1))

ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  mu[rand_idx,]<- runif(p)
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk)
  lld[i] <- mu_k.post(theta, mu, Sigma, Z, rand_idx)
}

lld-ldu

## Comparing Gibbs Sample and Posterior Distribution 
mu_o <- array(runif(K*p)*sample(1:5, K*p, replace=TRUE), dim =c(K, p))
init_vals<-list(rho=rho, V=V, Z=Z, theta=theta, Sigma=Sigma, mu = mu_o, vk = vk )
system.time({
  Theta<-gibbsSamples(y, 2e5, init_vals, burnin_period= 2e4, Sigma.fixed = TRUE, theta.fixed= TRUE, rho.fixed =TRUE, Z.fixed = TRUE)
})

hist(Theta$mu[,5,1], breaks = 100, freq = FALSE,
     xlab = expression(mu),
     main = expression(p(mu*" | "*Sigma,V,theta,rho,bold(z),bold(y))))
analytic_params <- mu_k_marginals.param(theta, 1, Sigma, Z, 5)
curve(dnorm(x = x, mean = analytic_params$mean, sd = analytic_params$sd), add = TRUE, col = "red")

hist(Theta$mu[,5,2], breaks = 100, freq = FALSE,
     xlab = expression(mu),
     main = expression(p(mu*" | "*Sigma,V,theta,rho,bold(z),bold(y))))
analytic_params <- mu_k_marginals.param(theta, 2, Sigma, Z, 5)
curve(dnorm(x = x, mean = analytic_params$mean, sd = analytic_params$sd), add = TRUE, col = "red")

hist(Theta$mu[,5,3], breaks = 100, freq = FALSE,
     xlab = expression(mu),
     main = expression(p(mu*" | "*Sigma,V,theta,rho,bold(z),bold(y))))
analytic_params <- mu_k_marginals.param(theta, 3, Sigma, Z, 5)
curve(dnorm(x = x, mean = analytic_params$mean, sd = analytic_params$sd), add = TRUE, col = "red")

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

rand_idx<-as.integer(N*runif(1))

ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  z_i <- Z[rand_idx]
  theta[,rand_idx]<- rmNorm(1, mu[z_i,], Sigma[,,z_i]) 
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk)
  lld[i] <- theta_i.post(y, theta, V, mu, Sigma, Z, rand_idx)
}

lld-ldu

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

rand_idx<-as.integer(K*runif(1))

ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  new_Sigma <-matrix(runif(p*p)*2-1, nrow= p, ncol =p ) 
  new_Sigma <- t(new_Sigma) %*% new_Sigma
  Sigma[,,rand_idx] <- new_Sigma 
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk)
  lld[i] <- sigma_k.post(theta, mu, Sigma, Z, rand_idx, vk)
}

lld-ldu

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

ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  rho[,1]<-runif(K)
  rho[,1]<-rho[,1]/sum(rho[,1])
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk)
  lld[i] <- rho.post(Z, rho[,1])
}

lld-ldu

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

rand_idx<-as.integer(N*runif(1))
print(rand_idx)
ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  Z[rand_idx]<- rcategorical(rho)
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk)
  lld[i] <- z.post(Z, rho[,1], theta, mu, Sigma, rand_idx)
}

lld-ldu

###################################  Testing Sampler ################################################

