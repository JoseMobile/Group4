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

fit1 <- gibbsSampler(data=theta_hat, V=obs_info, burnin_period=1e3, numIter=1e4, K=6, Theta.out=FALSE, z.out=FALSE)

fitted_clusters <- numeric(length=N)
for (i in 1:N){
  fitted_clusters[i] <- which.max(fit$Lambda[i, ])
}

err <- sum(cluster != fitted_clusters) / N

plot(cluster, fitted_clusters1)
