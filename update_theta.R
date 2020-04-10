require("mniw")

#' Conditional update of theta using the random effects multivariate normal
#' 
#' @param y A `n x p` matrix. Each row of y is an observation on p dimensions of the data
#' @param V A `p x p x n` array of covariance matricies for each observation of y.
#' @param mu A `K x p` matrix. The k-th row of mu is the mean of theta in class k. 
#' @param Sigma A `p x p x K` array of covariance matricies. The k-th matrix is the covariance matrix for theta in class k.
#' @param z A lenght `n` vector for the class observations of y.
#' @details `n` is the number of observations of y. `p` is the dimension of each theta, i.e. the number features of the data. `K` is the number of possible classes of y. 
#'
#'
update_theta <- function(y, V, mu, Sigma, z) {
  # dimension variables
  K <- nrow(mu) # classes
  n <- nrow(y)  # observations
  p <- ncol(y)  # features
  
  # create new mu and Sigma objects that match 
  # with y over the class observations
  
  # allocate space for new mu and Sigma
  new_mu <- matrix(nrow = n, ncol = p)
  new_Sigma <- array(dim = c(p,p,n))
  
  for (k in 1:K) {
    # iterate through classes and fill in new mu and Sigma with
    # rows of mu and Sigma based on observations from z
    
    ind <- which(z == k) # indices for z that are equal to k
    count <- length(ind) # number of observations of class k
    # fill in all of new_mu with class mean k with index in ind
    new_mu[ind,] <- matrix(mu[k,], nrow = count, ncol = p, byrow = TRUE)
    # fill in indicies of ind of new_Sigma with covariance matricies for k
    new_Sigma[,,ind] <- Sigma[,,k]
  }
  
  # sample using efficient random effects sampler with new mu and Sigma
  theta <- rRxNorm(n, y, V, new_mu, new_Sigma)
  return(theta)
}

# test function using parameters of the correct dimensions
#K <- 2
#p <- 4
#n <- 20

#y <- matrix(rnorm(n*p), nrow = n)
#V <- array(rnorm(p*p*n), dim = c(p,p,n))
#for (ii in 1:n) {
#  V[,,ii] <- crossprod(v[,,ii])
#}
#mu <- matrix(rnorm(K*p), nrow = K)
#Sigma <- array(rnorm(p*p*K), dim = c(p,p,K))
#for (k in 1:K) {
#  Sigma[,,k] <- crossprod(Sigma[,,k])
#}
#z <- ceiling(K*runif(n))

#update_theta(y, V, mu, Sigma, z)
