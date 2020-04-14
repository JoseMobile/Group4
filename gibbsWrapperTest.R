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

########### Testing Gibbs Sampler ####################
source('gibbs-wrapper.R')

system.time({
  Theta<-gibbsSampler(data=y, V=V, burnin_period=1, numIter=10, K=5, Theta.out=TRUE)
})

print(rowMeans(Theta$rho))



