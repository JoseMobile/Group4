
require(mniw)
require(numDeriv)
source('nnm-functions.R')


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
Omega<-data$Omega 

rand_idx<- as.integer(K*runif(1))

ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  mu[rand_idx,]<- runif(p)
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk, Omega)
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
Omega<-data$Omega 

rand_idx<-as.integer(N*runif(1))

ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  z_i <- Z[rand_idx]
  theta[rand_idx,]<- rmNorm(1, mu[z_i,], Sigma[,,z_i]) 
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk, Omega)
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
Omega<-data$Omega 

rand_idx<-as.integer(K*runif(1))

ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  new_Sigma <-matrix(runif(p*p)*2-1, nrow= p, ncol =p ) 
  new_Sigma <- t(new_Sigma) %*% new_Sigma
  Sigma[,,rand_idx] <- new_Sigma 
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk, Omega) 
  lld[i] <- sigma_k.post(theta, mu, Sigma, Z, rand_idx, vk, Omega)
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
Omega<-data$Omega 


ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  rho[,1]<-runif(K)
  rho[,1]<-rho[,1]/sum(rho[,1])
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk, Omega)
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
Omega<-data$Omega 


rand_idx<-as.integer(N*runif(1))
print(rand_idx)
ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  Z[rand_idx]<- rcategorical(rho)
  ldu[i]<- nnm_post(mu, Sigma, rho[,1], y, V, theta, Z, vk, Omega)
  lld[i] <- z.post(Z, rho[,1], theta, mu, Sigma, rand_idx)
}

lld-ldu

###################################  Testing Sampler ################################################


#finMix_mod <- stan_model("finiteMixtures.stan")
source('nnm-functions.R')
K<-5
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
Omega<-data$Omega 


for(i in 1:K){
  cat(i, ":", sum(i==Z), "\n")
}

cat(dim(y))
cat(length(Z))

mu_init <- array(runif(K*p), dim =c(K, p))
rho_init <- matrix(runif(K), nrow=K, ncol=1)
rho_init[,1]<-rho_init[,1]/sum(rho_init[,1])
theta_init <- array(runif(N*p), dim =c(N, p))
sigma_init <- array(runif(K*p*p)*2-1, dim =c(p,p,K))
for(i in 1:K){
  sigma_init[,,i]<-t(sigma_init[,,i]) %*% sigma_init[,,i]
}

Z_init <- rep(0,N)
  
while(TRUE){
  for (i in 1:N){
    Z_init[i]<-rcategorical(rho_init)
  }
  if (proper_sample(Z_init, K)) break
}

for(i in 1:K){
  cat(i, ":", sum(i==Z_init), "\n")
}

param_init = list(Sigma = sigma_init, mu= mu_init, theta=theta_init, rho=rho_init[,1], Z=Z_init, vk = vk, V= V, Omega= Omega)
#param_list = list(param_init)

system.time({
  Theta<-gibbsSamples(y, 1e4, param_init, burnin_period= 1e3)
})

print(rho)
print(Theta$mu[1,,])
print(colMeans(Theta$rho))

Theta$theta_samples



