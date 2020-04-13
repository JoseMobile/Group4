source('nnm-functions.R')
source('gibbs-wrapper.R')

K<-5
N<-100
p<-3

########### Real Data Simulation ###############
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

########### Starting Data Simulation ###############

# mu_init <- array(runif(K*p), dim =c(K, p))
# rho_init <- matrix(runif(K), nrow=K, ncol=1)
# rho_init[,1]<-rho_init[,1]/sum(rho_init[,1])
# theta_init <- array(runif(N*p), dim =c(N, p))
# sigma_init <- array(runif(K*p*p)*2-1, dim =c(p,p,K))
# for(i in 1:K){
#   sigma_init[,,i]<-t(sigma_init[,,i]) %*% sigma_init[,,i]
# }
# Z_init <- rep(0,N)
# 
# while(TRUE){
#   for (i in 1:N){
#     Z_init[i]<-rcategorical(rho_init)
#   }
#   if (proper_sample(Z_init, K)) break
# }
# 
# for(i in 1:K){
#   cat(i, ":", sum(i==Z_init), "\n")
# }

########### Testing Gibbs Sampler ####################
source('gibbs-wrapper.R')
system.time({
  Theta<-gibbsSampler(data=y, V=V, burnin_period=1e3, numIter=1e4, K=5, Theta.out=TRUE)
})


