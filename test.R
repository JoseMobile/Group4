
require(rstan)
require(mniw)
source('nnm-functions.R')

###################################  Simulate Data ################################################
#Helper Function to Simulate the Data 

sim_data<-function(K, N, p){
  K<-K
  N<-N
  p<-p 
  Z <- rep(0, N)
  
  rho<-matrix(rep(0, K), nrow=K, ncol=1)
  rho[,1]<-runif(K)
  rho[,1]<-rho[,1]/sum(rho[,1])
  
  for (i in 1:N){
    Z[i]<-rcategorical(rho)
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
  
  data = list(rho, y , V, Z, theta, mu, Sigma)
  names(data)<-c('rho', 'y', 'V', 'Z', 'theta', 'mu', 'Sigma')
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

rand_idx<-0 

while(TRUE){
  idx<- as.integer(K*runif(1))
  if (sum(idx == Z) > 0){
    rand_idx<-idx
    break
  }
}

ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  mu[rand_idx,]<- runif(p)
  ldu[i]<- nnm_loglik(mu, Sigma, rho[,1], y, V, theta, Z)
  lld[i] <- mu_k.loglik(theta, mu, Sigma, Z, rand_idx)
}

lld-ldu

###################################  Testing Theta ################################################

data<- sim_data(K, N, p)
rho <- data$rho
y <- data$y
V <-data$V
Z <- data$Z
theta <- data$theta
mu <-data$mu
Sigma <-data$Sigma

rand_idx<-as.integer(N*runif(1))

ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  z_i <- Z[rand_idx]
  theta[,rand_idx]<- rmNorm(1, mu[z_i,], Sigma[,,z_i]) 
  ldu[i]<- nnm_loglik(mu, Sigma, rho[,1], y, V, theta, Z)
  lld[i] <- theta_i.loglik(y, theta, V, mu, Sigma, Z, rand_idx)
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

rand_idx<-0

while(TRUE){
  idx<- as.integer(K*runif(1))
  if (sum(idx == Z) > 0){
    rand_idx<-idx
    break
  }
}


ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  new_Sigma <-matrix(runif(p*p)*2-1, nrow= p, ncol =p ) 
  new_Sigma <- t(new_Sigma) %*% new_Sigma
  Sigma[,,rand_idx] <- new_Sigma 
  ldu[i]<- nnm_loglik(mu, Sigma, rho[,1], y, V, theta, Z)
  lld[i] <- sigma_k.loglik(theta, mu, Sigma, Z, rand_idx)
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

ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  rho[,1]<-runif(K)
  rho[,1]<-rho[,1]/sum(rho[,1])
  ldu[i]<- nnm_loglik(mu, Sigma, rho[,1], y, V, theta, Z)
  lld[i] <- rho.loglik(Z, rho[,1])
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

rand_idx<-as.integer(N*runif(1))
print(rand_idx)
ldu<-rep(0, num_sim)
lld<-rep(0, num_sim)

for ( i in 1: num_sim){
  Z[rand_idx]<- rcategorical(rho)
  ldu[i]<- nnm_loglik(mu, Sigma, rho[,1], y, V, theta, Z)
  lld[i] <- z.loglik(Z, rho[,1], theta, mu, Sigma, rand_idx)
}

lld-ldu

