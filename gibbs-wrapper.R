#' Performs Gibbs Sampling given initial values
#' 
#' @param data The collected data for which we want to cluster the observations
#' @param M A scalar value indicating the desired number of MCMC samples to be generated
#' @param K A scalar value indicating the desired number of clusters 
#' @param M A scalar value indicating the desired number of MCMC samples to be generated


#' @param theta An 'n x p' matrix to act as initial value theta_0, where row i is the parameters of observation i.
#' @param mu An 'K x p' matrix to act as initial value mu_0, where the row k is expected value of theta if theta is in class k.
#' @param rho An 'K x 1' vector of class membership probabilities to act as initial value rho_0
#' @param Sigma A'p x p x K' array of variance-covariance matrices for each class k, to act as initial values Sigma_k,0
#' @param z A 'n x 1' vector of class memberships for each observation to act as initial value Z_0
#' 
#' @return 
#' 
#' @details 

gibbsSampler <- function(data, K, M=1000, theta=NULL, mu=NULL, rho=NULL, Sigma=NULL, z=NULL){
  # Check if arguments are NULL, and create initial values 
  if (is.null(mu)){
  
  }
  ...
  if (is.null(K)){
    ...
  }  
  
  N <- nrow(theta)
  total_Lambda <- matrix(0, nrow=N, ncol=K)
  
  # Create Defaults 
  
  # Need to choose how we want to do feature reduction 
  
  old_theta <- new_theta
  q <- ncol(old_theta)
  
  old_mu <- new_mu
  
  old_Sigma <- new_Sigma # Initialize some Sigma
  #  Maybe sigma <- matrix(rnorm(q * q), ncol=q, nrow=q), sigma <- t(sigma) %*% sigma
  
  old_rho <- rep(1/K, K)
  old_z <- sample(1:K, replace=TRUE) 
  
  # Compute V
  
  for (m in 1:M){
     # Call the update functions
     #  Check the dimensions of the return values match the other functions 
     
     # Update theta, 
     
     # Update mu, a 'K x p' matrix of class averages of theta
     new_mu <-  
     
     # Update Sigma
     
     # Update rho, a K x 1 vector
     new_rho <- update_rho(old_theta, old_z, old_mu, old_Sigma)
     
     # Update z, an N x 1 vector
     new_z <- update_z(old_theta, old_mu, old_rho, old_Sigma, give.Lambda=TRUE)
     
     # Update all values
     old_theta <- new_theta
     old_mu <- new_mu
     old_Sigma <- new_Sigma
     old_rho <- new_rho
     old_z <- new_z 
     
     total_Lambda <- total_Lambda + new_z$Lambda
  }
  
  total_Lambda <- total_Lambda / M
  
  # Compute class membership posterior probabilities 
  
  # Should we return the class memberships or a function that lets us do classification based on this "training set"
}
