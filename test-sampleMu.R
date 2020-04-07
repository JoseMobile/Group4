require("mniw")
source("sampleMu.R")
source("nnm-functions.R")

#' Function to simulate data for the Model-Based Clustering problem
#' 
#' @param K Number of classes. Scalar.
#' @param N Number of observations to simulate. Scalar.
#' @param p Number of features to include in simulation. Scalar.
sim_data<-function(K, N, p){
  K <- K
  N <- N
  p <- p 
  Z <- rep(0, N) # allocating space for Z vector

  # rho is the matrix of probablilities for the multinomial draws of Z
  # only a single column so that it works with the rcategorical function
  rho <- matrix(rep(0, K), nrow=K, ncol=1) # allocating space
  rho[,1] <- runif(K)             # filling in rho with random numbers
  rho[,1] <- rho[,1]/sum(rho[,1]) # normalize to probabilites that sum to one
  
  for (i in 1:N) { # fill in Z with multinomial draws
    Z[i] <- rcategorical(rho)
  }
  
  # allocate space for Sigma and fill with random numbers
  Sigma <- array(runif(p*p*K)*2-1, dim=c(p,p,K)) # pxpxK array of numbers in (1,2)
  for (i in 1:K) { # iterate through classes to make Sigma symmetric
    Sigma[,,i]<- t(Sigma[,,i]) %*% Sigma[,,i] # Matrix %*% Matrix' is symmetric
  }
  
  # allocate space for V and fill with random numbers
  V <- array(runif(p*p*N)*2-1, dim=c(p,p,N)) # pxpxN array of numbers in (1,2)
  for (i in 1:N) { # iterate through sample size to make V symmetric 
    V[,,i]<- t(V[,,i]) %*% V[,,i] # V %*% V' is symmetric
  }
  
  # allocate space for mu and theta matricies
  mu <- matrix(runif(K*p), nrow = K, ncol = p)   # mu is Kxp, random numbers in (0,1)
  theta <- matrix(runif(0, p*N), nrow =p, ncol=N) # theta is pXN, filled with zeros
  
  for(i in 1:N) { # filling in theta matrix with matrix normal
    z_i <- Z[i]   # distribution is based on class observation
    theta[,i] <- rmNorm(n=1, mu[z_i,], Sigma[,,z_i]) # sample
  }
  
  # sampling for y vectors
  y<- matrix(rep(0, p*N), nrow =N, ncol=p) # allocating space for Nxp matrix
  for(i in 1:N) { 
    y[i,]<- rmNorm(n=1, theta[,i], V[,,i])
  }
  
  # create data object to return all values
  data = list(rho, y , V, Z, theta, mu, Sigma)
  names(data)<-c('rho', 'y', 'V', 'Z', 'theta', 'mu', 'Sigma')
  data # return data object
}


# See if sample Mu is close to mean of theta ----------------------------------


K <- 3 # classes
q <- 10 # features
n <- 1000 # observations

# simulate data and set to easier variable names
simData <- sim_data(K, n, q)
theta <- t(simData$theta) # sim_data gives 'qxn' where sample mu uses 'nxq'
Sigma <- simData$Sigma
z <- simData$Z

# observed theta that mu is supposed to me mean of
meanTheta.hat <- matrix(nrow = K, ncol = q)
for (k in 1:K) { # iterate thorugh classes
  meanTheta.hat[k,] <- colMeans(subset(theta, z == k))
}

# going to test if sample mu is sampling reasonable means by using a 
# hypothesis test on each variable after uncorrelating the multivariate
# normal distibution
num_sim <- 1 # times to repeat test

# observed theta that mu is supposed to me mean of
theta_mean <- matrix(nrow = K, ncol = q)
for (k in 1:K) { # iterate thorugh classes
  theta_mean[k,] <- colMeans(subset(theta, z == k))
}

# inverse of cholesky decomposition of Sigma
chol_inv <- Sigma # allocating space for chol_inv
for (k in 1:K) {
  chol_inv[,,k] <- chol(Sigma[,,k])
  chol_inv[,,k] <- solve(chol_inv[,,k])
}

# hypothesis testing if mu could be mean of thetas
hyp_test <- array(dim = c(K, q, num_sim)) # array to hold hypothesis tests
save_par <- par()
par(mfrow = c(K,q/2)) # set grid of graphs so can show all histograms
for (rep in 1:num_sim) { # find hypothesis test values
  mu <- sampleMu(theta, Sigma, z, currMu = matrix(0, K, q))
  
  for (k in 1:K) { # iterate through rows of mu, classes
    mu_trans <- chol_inv[,,k] %*% (simData$mu[k,] - theta_mean[k,]) # transform mu to be uncorrelated
    # mu_trans should be distributed iid N(0,1) if the code is correct
    # hypothesis test if mu is close to zero
    hyp_test[k,,rep] <- 2*dnorm(abs(mu_trans)) # pvalue for N(0,1)
    
    for (jj in 1:q) { # create histograms to see if mu could be mean
      # seeing if sampled mu is reasonable for mean of theta
      #hist(subset(theta[,jj], z == k))
      #abline(v = mu[k,jj], col = "red")
    }
  }
}
par(save_par) # reset graphical parameters

# test how many means are rejected by the hypothesis test
rjt_test <- hyp_test <= 0.05 # rejected at 0.05 level
sum(rjt_test)/length(rjt_test) # percentage rejected

################################################################################
# Testing if there are no observations of a class ---------------------------------
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

