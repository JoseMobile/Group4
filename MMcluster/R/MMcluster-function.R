#' Mixture Model Clustering
#'
#' @description Performs mixture model clustering using a Gibbs sampler by assuming the data follows a random effects mixture model.
#'
#' @param data Input matrix with dimensions `N x p`.
#' @param V The observed Fisher information matrix of each observation. Given as array with dimensions `p x p x N`.
#' @param prior An optional list containing the parameters for the prior selected. Must have elements named `Vk` and `Omega`.
#'              `Vk` being a vector of length `K` and `Omega` being an array of variance matricies with dimension `p x p x K`.
#' @param initParamVals A list containing the initial values of the parameters of the model. Must have elements names `theta`, `mu`, `Sigma`, `rho`, and `z`. Where
#'                      `theta` is a `N x p` matrix,
#'                      `mu` is a `K x p` matrix,
#'                      `Sigma` is a `p x p x K` array of variance matricies,
#'                      `rho` is a probability vector with length `K`,
#'                      and `z` is a vector of length `N` with each element being an integer within \{1,...,K\}.
#' @param K A scalar value indicating the desired number of clusters.
#' @param burnin_period A scalar value indicating the initial number of samples to discard.
#' @param numIter A scalar value indicating the desired number of MCMC samples to be generated.
#' @param mu.fixed A true or false indicating whether to sample new mu values. See details for the usage of the mu parameter.
#' @param Sigma.fixed  A true or false indicating whether to sample new Sigma values. See details for the usage of the Sigma parameter.
#' @param theta.fixed A true or false indicating whether to sample new theta values. See details for the usage of the theta parameter.
#' @param rho.fixed A true or false indicating whether to sample new rho values. See details for the usage of the rho parameter.
#' @param z.fixed A true or false indicating whether to sample new z values. See details for the usage of the z parameter.
#' @param Theta.out A true or false indicating whether to save and return all Theta samples. Caution: This will be a very large object.
#' @param z.out A true or false indicating whether to save all z samples. Caution: This will be a very large object.
#' @return A list of all the sampled parameters to the posterior distribution of z. The named elements are as follows: \cr
#'         `theta` is a `N x p` matrix or Theta.out = TRUE a `N x p x M` array. \cr
#'         `mu` is a `K x p x M` array. \cr
#'         `Sigma` is a `p x p x K x M` array of variance matricies. \cr
#'         `rho` is a `K x M` matrix with the rows being probability vectors (sum to 1). \cr
#'         `z` is a length `N` vector with each element being an integer within \{1,...,K\} or if z.out = TRUE a `N x M` matrix.
#' @details The hierarchical model used in the Gibbs sampler is as follows:
#'          \deqn{z ~ Multinom(K, \rho)}
#'          Where K is the number of clusters and \eqn{\rho} is a probability vector for being in each group. Each z is identically independently distributed.
#'          \deqn{\theta | z ~ N(\mu, \Sigma)}
#'          Where there are K different \eqn{\mu} and \eqn{\Sigma} parameters, one each for every cluster. This cluster is given by the cluster observation of z.
#'          \deqn{y | \theta ~ N(\theta, V)}
#'          Where y is the data given to the function and V is the Fisher information matrix of each observation. So, \eqn{\theta} creates a random effects model of the data for y.
#' @export
MMcluster <- function(data, V, prior, initParamVals, K, burnin_period, numIter, mu.fixed = FALSE,
                      Sigma.fixed =FALSE, theta.fixed= FALSE, rho.fixed = FALSE, z.fixed = FALSE, Theta.out=FALSE, z.out=FALSE){
  N <- nrow(data)
  p <- ncol(data)

  #Getting the Prior Parameters
  if (missing(prior)){
    Vk <- rep(p+2,K)
    Omega<- array(rep(0, p*p*K), dim=c(p,p,K))

    for(i in 1:K)
      Omega[,,i]<- var(data) *(N-1)/(N-p)
  } else {
    if ( (! "Vk" %in% names(prior) ) || (!"Omega" %in% names(prior)) )
      stop("prior needs Vk and Omega params!")

    Vk<- prior$Vk
    Omega<- prior$Omega
  }

  #Getting the Initial Value of Model Parameters
  if (missing(initParamVals)){ # no initial values given
    # INITALIZING PARAMETERS

    # Compute overall data means and variances
    data_mu <- colMeans(data)
    data_var <- var(data) * (N-1)/(N-p) # variance of data

    # initalizing parameters using kmeans++ algorithm
    old_z <- init_z(data, K) # use kmeans++ algorith to initialize z

    # set rho to proportion of each observation
    count <- as.numeric(table(old_z))

    # table gives count of each observation of class
    # then dividing by N gives the proportion
    old_rho <- count/N

    # set theta to observed data
    old_theta <- data

    # allocate space for mu and Sigma
    old_mu <- matrix(nrow = K, ncol = p)
    old_Sigma <- array(dim = c(p,p,K))
    # set mu to sample mean of data in that class
    for (k in 1:K) {
      ind <- which(old_z == k) # indices of data in class k
      if (count[k] > p) { # enough data to sample both mu and Sigma
        old_mu[k,] <- colMeans(data[ind,]) # sample mu
        # sample Sigma, var() scaled by n-1, so rescale to n-p
        old_Sigma[,,k] <- var(data[ind,]) * (count[k]-1)/(count[k]-p)
      } else {
        # not enough data for Sigma, so use data variance
        old_Sigma[,,k] <- data_var
        # sampling for mu
        if (count[k] > 1) { # enough data for sampling mu
          old_mu[k,] <- colMeans(data[ind,])
        } else if (count[k] == 1) { # only one so set to mu
          old_mu[k,] <- data[ind,]
        } else { # count is zero, so use data mean
          old_mu[k,] <- data_mu
        } # end second if
      } # end first if
    } # end for

  } else { # given inital parameters
    if ( (! "theta" %in% names(initParamVals) ) || (!"mu" %in% names(initParamVals)) ||  (!"Sigma" %in% names(initParamVals))
         ||  (!"rho" %in% names(initParamVals))||  (!"z" %in% names(initParamVals)))
      stop("initParamVals needs theta, mu, Sigma, rho and z params!")

    old_theta<- initParamVals$theta
    old_mu<- initParamVals$mu
    old_Sigma<- initParamVals$Sigma
    old_rho<- initParamVals$rho
    old_z<- initParamVals$z
  }

  # Allocate space to store all generated samples
  rho <- array(NA, dim=c(K, numIter))
  mu <- array(NA, dim=c(K, p, numIter))
  Sigma <- array(NA, dim=c(p, p, K, numIter))
  if (Theta.out)
    Theta <- array(NA, dim=c(p, N, numIter))
  else
    Theta<- NULL

  if (z.out)
    z <- array(NA, dim=c(N, numIter))
  else
    z<- NULL

  total_Lambda <- array(0, dim=c(N, K))

  # The dimensions of the arguments are:
  #   theta: N x p      mu: K x p     Sigma: p x p x K    rho: K x 1    z: N x 1
  # Parameters orders should all be: theta, mu, Sigma, rho, z
  #   If not check the function calls and standardize them
  for (m in 1:(burnin_period +numIter)){
    # output iteration number every 200 runs
    if (m %% 200 == 0) {
      cat("Iteration # ", m, "\n")
    }

    # Call the update functions

    if (!theta.fixed){
    # Update theta, a 'N x p' matrix of random effects of y
      new_theta <- update_theta(data, V, old_mu, old_Sigma, old_z)
    }
    if (!mu.fixed){
    # Update mu, a 'K x p' matrix of class averages of theta
      new_mu <- update_mu(new_theta, old_mu, old_Sigma, old_z)
    }
    if(!Sigma.fixed){
    # Update Sigma, a 'p x p x K' array of covariance matricies
      new_Sigma <- update_Sigma(data, new_theta, new_mu, old_Sigma, old_z, Omega, Vk)
    }
    if(!rho.fixed){
    # Update rho, a K x 1 vector
      new_rho <- update_rho(new_theta, new_mu, new_Sigma, old_z)
    }
    if (!z.fixed){
    # Update z, an N x 1 vector
      new_z <- update_z(new_theta, new_mu, new_Sigma, new_rho)
    }

    # Store theta and z values if desired
    if (m > burnin_period){
      curIdx<- m - burnin_period
      if (Theta.out){
        Theta[, , curIdx] <- new_theta
      }
      if (z.out){
        z[, curIdx] <- new_z$z
      }

      rho[ ,curIdx]<- new_rho
      mu[ , ,curIdx]<- new_mu
      Sigma[ , , ,curIdx]<-  new_Sigma
      total_Lambda <- total_Lambda + new_z$Lambda
    }

    # Update all values
    old_theta <- new_theta
    old_mu <- new_mu
    old_Sigma <- new_Sigma
    old_rho <- new_rho
    old_z <- new_z$z
  }

  total_Lambda <- total_Lambda / numIter

  return(list(theta=Theta, mu=mu, Sigma=Sigma, rho=rho, z=z, Lambda=total_Lambda))
}
