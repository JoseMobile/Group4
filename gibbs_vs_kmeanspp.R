source("gibbs-wrapper.R")
require("plot3D")
# install.packages("plot3D") # sometimes require() does not install the plot3D package, if so, run this line

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

# getting true cluster affiliations
real_N_k <- numeric(length=K)
for (k in 1:K){
  real_N_k[k] <- sum(cluster == k)
}

# fit gibbs sampler
gibbs_fit <- gibbsSampler(data=theta_hat, V=obs_info, burnin_period=1e3, numIter=1e4, K=6, Theta.out=FALSE, z.out=FALSE)

# get estimated cluster affiliation from gibbs sampler
fitted_clusters <- numeric(length=N)
for (i in 1:N){
  fitted_clusters[i] <- which.max(gibbs_fit$Lambda[i, ])
}

# error rate
err <- sum(cluster != fitted_clusters) / N
err

# set colors for plotting
index <- 1:N
gibb_cols <- (cluster != fitted_clusters)
gibb_cols[gibb_cols == TRUE] <- "red"
gibb_cols[gibb_cols == FALSE] <- "green"

# plots to visualize results
#plot(cluster, fitted_clusters,col = gibb_cols, xlab = "true_cluster",ylab="fitted_clusters",main ="real cluster vs fitted")
scatter3D(x = cluster,y=fitted_clusters,z = index,col = gibb_cols,xlab="true_cluster",ylab="fitted_cluster",zlab="index",phi=-20,theta=-30,main = "real cluster vs fitted")

# clustering real data with kmeans++ algorithm
# note since kmeans starts with random centroids, individual kmeans can have different results.
# We will conduct many kmeans and aggregate the results

#' Initial cluster allocation.
#'
#' Initializes the clusters using the kmeans++ algorithm of Arthur & Vassilvitskii (2007).
#'
#' @param y An `N x p` matrix of observations, each of which is a column.
#' @param K Number of clusters (positive integer).
#' @param max_iter The maximum number of steps in the [stats::kmeans()] algorithm.
#' @param nstart The number of random starts in the [stats::kmeans()] algorithm.
#'
#' @return A vector of length `N` consisting of integers between `1` and `K` specifying an initial clustering of the observations.
#' @author Martin Lysy
#' @references Arthur, D., Vassilvitskii, S. "k-means++: the advantages of careful seeding" *Proceedings of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms.* (2007): 1027–1035. <http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf>.
init_z <- function(y, K, max_iter = 10, nstart = 10) {
  # init++
  N <- nrow(y)
  p <- ncol(y)
  x <- t(y) # easier to use columns
  centers <- matrix(NA, p, K) # centers
  icenters <- rep(NA, K-1) # indices of centers
  minD <- rep(Inf, N) # current minimum distance vector
  # initialize
  icenters[1] <- sample(N,1)
  centers[,1] <- x[,icenters[1]]
  for(ii in 1:(K-1)) {
    D <- colSums((x - centers[,ii])^2) # new distance
    minD <- pmin(minD, D)
    icenters[ii+1] <- sample(N, 1, prob = minD)
    centers[,ii+1] <- x[,icenters[ii+1]]
  }
  centers <- t(centers)
  colnames(centers) <- rownames(x)
  # run kmeans with these centers
  km <- kmeans(x = y, centers = centers, nstart = nstart, iter.max = max_iter)
  km$cluster
}

max_iters <- 10000  # set number of kmeans to run
kmeans_results <- matrix(nrow = N,ncol = max_iters) # pre-allocate matrix

# simple function to calculate the mode in vector of numbers
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# get results from max_iters kmeans
for (iter in 1:max_iters){
  kmeans_results[,iter] <- init_z(y = theta_hat,K = 6)
}

fitted_clusters1 <- matrix(nrow = N,ncol = 1) # pre-allocate matrix to store estimated clusters

# the mode cluster membership for each observation becomes the estimate cluster membership
for (row in 1:N){
  fitted_clusters1[row,] <- Mode(kmeans_results[row,])
}

err1 <- sum(cluster != fitted_clusters1) / N
err1
# kmeans++ x 10000 performs badly (> 75% error) but still performed better than the gibbs sampler (> 80% error) 
# setting colors for plotting
kmeans_cols <- (cluster != fitted_clusters1)
kmeans_cols[kmeans_cols == TRUE] <- "red"
kmeans_cols[kmeans_cols == FALSE] <- "green"

# plots to visualize results
#plot(cluster, fitted_clusters1,col = gibb_cols, xlab = "true_cluster",ylab="fitted_clusters",main ="real cluster vs fitted (kmeans++ x 10000)")
scatter3D(x = cluster,y=fitted_clusters1,z = index,col = kmeans_cols,xlab="true_cluster",ylab="fitted_cluster",zlab="index",phi=-20,theta=-30,main = "real cluster vs fitted (kmeans++ x 10000)")

#-----------------------------------------------------------------------------------------------------------
# testing with standardized data (data centered at 0 with sd of 1)
#-----------------------------------------------------------------------------------------------------------
theta_hat_scaled <- scale(theta_hat)

# fit gibbs sampler
gibbs_fit1 <- gibbsSampler(data=theta_hat_scaled, V=obs_info, burnin_period=1e3, numIter=1e4, K=6, Theta.out=FALSE, z.out=FALSE)

# get estimated cluster affiliation from gibbs sampler
fitted_clusters2 <- numeric(length=N)
for (i in 1:N){
  fitted_clusters2[i] <- which.max(gibbs_fit1$Lambda[i, ])
}

# error rate
err2 <- sum(cluster != fitted_clusters2) / N
err2

# set colors for plotting
gibb_cols1 <- (cluster != fitted_clusters2)
gibb_cols1[gibb_cols1 == TRUE] <- "red"
gibb_cols1[gibb_cols1 == FALSE] <- "green"

# plots to visualize results
#plot(cluster, fitted_clusters,col = gibb_cols, xlab = "true_cluster",ylab="fitted_clusters",main ="real cluster vs fitted")
scatter3D(x = cluster,y=fitted_clusters2,z = index,col = gibb_cols1,xlab="true_cluster",ylab="fitted_cluster",zlab="index",phi=-20,theta=-30,main = "real cluster vs fitted")

max_iters <- 10000
kmeans_results1 <- matrix(nrow = N,ncol = max_iters)

# get results from max_iters kmeans
for (iter in 1:max_iters){
  kmeans_results1[,iter] <- init_z(y = theta_hat_scaled,K = 6)
}

fitted_clusters3 <- matrix(nrow = N,ncol = 1) # pre-allocate matrix

# the mode cluster membership for each observation becomes the estimate cluster membership
for (row in 1:N){
  fitted_clusters3[row,] <- Mode(kmeans_results1[row,])
}

err3 <- sum(cluster != fitted_clusters3) / N
# kmeans++ x 10000 performs badly (> 75% error) but still performed similar to the gibbs sampler  
# setting colors for plotting
kmeans_cols1 <- (cluster != fitted_clusters3)
kmeans_cols1[kmeans_cols1 == TRUE] <- "red"
kmeans_cols1[kmeans_cols1 == FALSE] <- "green"

# plots to visualize results
#plot(cluster, fitted_clusters1,col = gibb_cols, xlab = "true_cluster",ylab="fitted_clusters",main ="real cluster vs fitted (kmeans++ x 10000)")
scatter3D(x = cluster,y=fitted_clusters3,z = index,col = kmeans_cols1,xlab="true_cluster",ylab="fitted_cluster",zlab="index",phi=-20,theta=-30,main = "real cluster vs fitted (kmeans++ x 10000)")

#--------------------------------------------
# clustering with simulated data

n <- 600
A <- runif(n = 200,min = -1, max = 1)
B <- runif(n = 200, min = 100, max = 110)
C <- runif(n = 200, min = -200, max = -190 )
D <- runif(n = n, min = -5, max = 5)

sim_data <- matrix(nrow = n,ncol = 2)
sim_data[,1] <- c(A,B,C)
sim_data[,2] <- D
var_sim1 <- var(sim_data[,1])
var_sim2 <- var(sim_data[,2])
cov_sim <- cov(sim_data[,1],sim_data[,2])
sim_obs_info <- solve(matrix(data = c(var_sim1,cov_sim,
                                cov_sim,var_sim2),nrow = 2,ncol=2,byrow=TRUE)) # inverse of fisher obvs is the asymptotic covariance matrix,
                                                                              # we will use the inverse of the observed covariance matrix as an estimate for the fisher obvs


plot(x=sim_data[,1],y=sim_data[,2])

# fit gibbs sampler
gibbs_fit2 <- gibbsSampler(data=sim_data, V=sim_obs_info, burnin_period=1e3, numIter=1e4, K=3, Theta.out=FALSE, z.out=FALSE)

# get estimated cluster affiliation from gibbs sampler
fitted_clusters4 <- numeric(length=n)
for (i in 1:600){
  fitted_clusters4[i] <- which.max(gibbs_fit2$Lambda[i, ])
}

# indices
index1 <- 1:n
# error rate
sim_results <- as.data.frame(sim_data)
sim_results$col <- fitted_clusters4
sim_results$col[sim_results$col == 1] <- "red"
sim_results$col[sim_results$col == 2] <- "green"
sim_results$col[sim_results$col == 3] <- "blue"

plot(x=sim_results$V1,y = sim_results$V2,col = sim_results$col, main="clusters determined by gibbs sampler")

