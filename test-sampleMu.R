require("mniw")
source("sampleMu.R")

# sample from a multi-variate normal using R code and see 
# if it is close to the answer from sampleMu

K <- 2  # number of classes
q <- 2  # number of features of theta
n <- 5 # sample size

# random vector for the classes of observations
z <- ceiling(K*runif(n))

# class means of each feature of theta
# rows are each class, each column is a feature
classMeanTheta <- matrix(3*runif(q*K), nrow = K, ncol = q)

# allocating space for Sigma
Sigma <- array(dim = c(q,q,K))

# fill in Sigma with random variance matricies
for (kk in 1:K) {
  Sigma[,,kk] <- rnorm(q*q)   # random numbers temporarily
  # multiply by transpose to make symmetric matrix which is random
  Sigma[,,kk] <- Sigma[,,kk] %*% t(Sigma[,,kk]) 
}

# create a random generated theta
theta <- matrix(nrow = n, ncol = q)
for (ii in 1:n) {
  k <- z[ii] # class of observation ii
  C <- chol(Sigma[,,k]) # cholesky decomp: C'C = Sigma
  Z <- rnorm(q) # random N(0,1) observations
  theta[ii,] <- Z %*% C + classMeanTheta[k,] # transform Z to N(meanTheta, C'C)
}

colMeans(subset(theta, z == 1))
colMeans(subset(theta, z == 2))

# test if sampleMu is close to generated classMeanTheta using the bootstrap
B <- 100 # number of bootstrap samples
muBootArray <- array(dim = c(q,q,B)) # allocating space

# bootstrap loop
for (b in 1:B) {
  muBootArray[,,b] <- sampleMu(theta, Sigma, z) # generate each sample
}

muBoot <- matrix(nrow = q, ncol = q) # allocating space for mean of bootstrap
muBootSD <- matrix(nrow = q, ncol = q) # allocating space for sd of bootstrap
for (ii in 1:q) {   # iterate through rows
  for (jj in 1:q) { # iterate through columns
    muBoot[ii,jj] <- mean(muBootArray[ii,jj,]) # mean of bootstrap
    muBootSD[ii,jj] <- sd(muBootArray[ii,jj,]) # sd of bootstrap
  }
}

# observed theta that mu is based off of
meanTheta.hat <- matrix(nrow = K, ncol = q)
for (k in 1:K) { # iterate thorugh classes
  meanTheta.hat[k,] <- colMeans(subset(theta, z == k))
}

# test to see if mean of theta is within mu bootstrap 95% CI
muBoot.lower <- muBoot - meanTheta.hat - 1.96 * muBootSD
muBoot.upper <- muBoot - meanTheta.hat + 1.96 * muBootSD


# Testing if there are no observations of a class
# function should return the same row for that class as was given in currMu

K <- 3  # number of classes
q <- 2  # number of features of theta
n <- 20 # sample size

# random vector for the classes of observations
z <- ceiling((K-1)*runif(n)) # The Kth class will never show up

# class means of each feature of theta
# rows are each class, each column is a feature
classMeanTheta <- matrix(3*runif(q*K), nrow = K, ncol = q)

# allocating space for Sigma
Sigma <- array(dim = c(q,q,K))

# fill in Sigma with random variance matricies
for (kk in 1:K) {
  Sigma[,,kk] <- rnorm(q*q)   # random numbers temporarily
  # multiply by transpose to make symmetric matrix which is random
  Sigma[,,kk] <- Sigma[,,kk] %*% t(Sigma[,,kk]) 
}

# create a random generated theta
theta <- matrix(nrow = n, ncol = q)
for (ii in 1:n) {
  k <- z[ii] # class of observation ii
  C <- chol(Sigma[,,k]) # cholesky decomp: C'C = Sigma
  Z <- rnorm(q) # random N(0,1) observations
  theta[ii,] <- Z %*% C + classMeanTheta[k,] # transform Z to N(meanTheta, C'C)
}

# sample mu by giving a matrix of zeros as the current mu
# this should have the last row (K) should always be zero
currMu <- matrix(rnorm(K*q), K, q) # Kxq matrix to give to sampleMu
Mu <- sampleMu(theta, Sigma, z, currMu)
Mu[K,] - currMu[K,] # should be zero