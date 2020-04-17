#' Conditional update of theta using the random effects multivariate normal
#'
#' @param y A `n x p` matrix. Each row of y is an observation on p dimensions of the data
#' @param V A `p x p x n` array of covariance matricies for each observation of y.
#' @param mu A `K x p` matrix. The k-th row of mu is the mean of theta in class k.
#' @param Sigma A `p x p x K` array of covariance matricies. The k-th matrix is the covariance matrix for theta in class k.
#' @param z A lenght `n` vector for the class observations of y.
#' @details `n` is the number of observations of y. `p` is the dimension of each theta, i.e. the number features of the data. `K` is the number of possible classes of y.
#' @return A `n x p` matrix. A random sample from the random effects multivariate normal distribution. If the return value is called theta then the details are as follows:
#' \deqn{\theta ~ N(\mu, \Sigma)}
#' \deqn{y | \theta ~ N(\theta, V)}
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
  count <- numeric(length = K)


  for (k in 1:K) {
    # iterate through classes and fill in new mu and Sigma with
    # rows of mu and Sigma based on observations from z

    ind <- which(z == k) # indices for z that are equal to k
    count[k] <- length(ind) # number of observations of class k
    # fill in all of new_mu with class mean k with index in ind
    # only execute if there are class observations to avoid warnings
    if (count[k] > 0) {
      new_mu[ind,] <- matrix(mu[k,], nrow = count[k], ncol = p, byrow = TRUE)
      # fill in indicies of ind of new_Sigma with covariance matricies for k
      new_Sigma[,,ind] <- Sigma[,,k]
    }
  }

  # sample using efficient random effects sampler with new mu and Sigma
  theta <- mniw::rRxNorm(n, y, V, new_mu, new_Sigma)
  return(theta)
}



#' Conditional update for `mu_k` by randomly sampling from the Multivariate Normal model with mean theta and variance proportional to Sigma (see details).
#'
#'
#' @param theta Estimated random effects of the mixture-normal model. A `n x p` matrix.
#' @param old_mu Current estimated mu. A `K x p` matrix.
#' @param Sigma Variance matricies of theta. A `p x p x K` array.
#' @param z Estimated classes of the y values. A `n x 1` vector.
#' @details `K` is the number of classes. `p` is the dimension of one observations of `theta_i`, i.e. the number of features. `n` is the number of observations. If a class `k` does not appear in the z vector, the `k`th row of the return matrix is the `k`th row of currMu. If a class has no observations in the z vector this mu is not updated because having a count of zero causes a division by zero.
#' @return A `K x 1` vector of randomly sampled `mu_k` values from the distribution `mu_k ~ Normal(theta_k, Sigma_k/N_k)`
#' ```
#' where `theta_k` is the mean of all `theta_i`s where `z_i = k` and `N_k` is the number of `z_i = k` (sum of indicator `z_k = k`).
update_mu <- function(theta, old_mu, Sigma, z) {
  # number of classes, K. number of features, q
  SigmaDim <- dim(Sigma)
  K <- SigmaDim[3]  # q x q x K
  q <- SigmaDim[1]  # q x q x K

  # allocating space for class mean theta to be used
  classMeanTheta <- matrix(nrow = K, ncol = q) # Means, theta_k

  # if there are no observations of a class k in z, the kth row of Mu won't change
  # after the for loop Mu[k,] is set to currMu[k,]
  # use boolean variable to simplify the process
  replaceMu <- rep(FALSE, K) # initalize to false, changed in for loop if needed

  # iterate through the classes
  for (kk in 1:K) {
    zInd <- (z == kk)  # indicator variable if z_i = k
    # have to check how many observations of class in z
    count <- sum(zInd)

    # Class mean theta has rows of classMeanTheta filled with mean of class
    if (count > 1) { # more than one observation of class
      classMeanTheta[kk,] <- colMeans(theta[zInd,])

      # Calculating scaled Sigma, divided by count of kk in z
      Sigma[,,kk] <- Sigma[,,kk]/count
    } else if (count == 1) { # only one observation of class
      classMeanTheta[kk,] <- theta[zInd,]

      # Calculating scaled Sigma, divided by count of kk in z
      Sigma[,,kk] <- Sigma[,,kk]/count
    } else { # no observations of class
      classMeanTheta[kk,] <- rep(0, q) # fill in theta so that can sample
      replaceMu[kk] <- TRUE # boolean variable used after samplin to replace row
    }


  }

  Mu <- mniw::rmNorm(K, classMeanTheta, Sigma) # sample Mu

  # replace rows with no observations in that class with last row mean
  Mu[replaceMu,] <- old_mu[replaceMu,]

  return(Mu)
}


#' Samples an array of matrices denoting the covariance matrices for each cluster
#'
#'
#' @param y A `n x p` matrix. Each row of y is an observation on p dimensions of the data.
#' @param theta A `n x p` matrix where each row is an observation and the columns are a subset of a larger parameter set.
#' @param mu An `K x p` vector denoting the prior cluster means.
#' @param old_Sigma A `p x p x K` array of variance matricies. The previous observation of Sigma.
#' @param z A length `n` vector denoting the class that each row of theta belongs to.
#' @param Omega The covariance matrixes of our chosen prior.
#' @param vK The virtual counts corresponding to mu of chosen prior
#'
#' @return A `p x p x K` array of matrices sampled from the inverse-wishart distribution where the kth matrix
#' denotes the posterior covariance matrix of the kth cluster.
#' @details The array of matrices that is returned is  `Sigma_k | A \ \{theta_k\}` which had been sampled from
#' ```
#' `Inverse-Wishart(\sum{\Theta_i - mu_k},N + vK)` where N is a vector containing cluster counts
#' , p is a vector of the sizes of theta and 1 is a vector of 1's. Assume that the dimension requirements
#' are met
#' ```
update_Sigma <- function(y, theta, mu, old_Sigma, z, Omega, vK) {
  # iterate through classes and make a new mu to match theta
  dim_mu <- dim(mu) # dimensions of mu c(K, p)
  K <- dim_mu[1] # rows are K
  p <- dim_mu[2] # columns are p
  n <- nrow(theta) # rows are n, number of observations

  # allocate space before
  replaceSigma <- rep(FALSE, K) # indicator for if there are empty clusters
  # initialized at false, will change if empty cluster if found
  count <- rep(NA, K) # allocate space for count of classes in z
  Psi <- array(dim = c(p,p,K)) # allocate space for Psi matrix

  for (k in 1:K) {
    k_inds <- which(z == k) # indices of observations of class k
    # calculate count, number of class observations in z
    count[k] <- length(k_inds)

    # initiate Psi_k to Omega_k
    Psi[,,k] <- Omega[,,k]

    # check how many observations there are in class k
    if (count[k] == 0) {
      replaceSigma[k] <- TRUE # found empty cluster
    } else { # cluster has at least one observation
      for (ii in 1:count[k]) { # iterate through all of the class observations
        ind <- k_inds[ii] # ii-th observation of class k

        # calculate residual of theta and mu
        res <- theta[ind,] - mu[k,]
        Psi[,,k] <- Psi[,,k] + tcrossprod(res) # outer product
        # tcrossprod(x) is x %*% t(x)
      }
      # Psi_k is now the sum of outer products of the rows of res
    }
  }

  # calculate nu
  nu <- count + vK # a vector of length K
  Sigma <- mniw::riwish(K, Psi, nu)

  # replace empty cluster class Sigmas with the previous Sigma
  Sigma[,,replaceSigma] <- old_Sigma[,,replaceSigma]

  return(Sigma)
}

#' Samples cluster memberships probabilities from a Dirichlet distribution
#'
#'
#' @param theta An 'n x p' matrix, where row i is the parameters of observation i.
#' @param mu An 'K x p' matrix where the row k is expected value of theta if theta is in class k.
#' @param Sigma A 'p x p x K' array of variance-covariance matrices for each class k.
#' @param Z An 'n x 1' vector of class memberships for each observation
#'
#' @return An 'n x 1' vector of updated class membership probabilities.
update_rho <- function(theta, mu, Sigma, Z){
  # Initialize constants and pre-allocate memory
  K <- nrow(mu)
  N <- nrow(theta)

  # Create alpha vector of class counts
  alpha <- numeric(length=K)
  for (k in 1:K){
    alpha[k] <- sum(Z == k) + 1
  }

  rho<-rdirichlet(1, alpha)
  return(rho)
}


#' Samples group memberships from a multinomial distribution conditioned on the other parameters
#'
#'
#' @param theta An 'n x p' matrix, where row i is the parameters of observation i.
#' @param mu An 'K x p' matrix where the row k is expected value of theta if theta is in class k.
#' @param Sigma A'p x p x K' array of variance-covariance matrices for each class k.
#' @param rho An 'K x 1' vector of class membership probabilities
#'
#' @return An 'n x 1' vector of updated class memberships, or if give.Lambda=TRUE, returns a list that contains the 'n x 1' updated class memberships and the computed 'n x k' Lambda matrix.
#'
update_z <- function(theta, mu, Sigma, rho){
  # Initialize constants and pre-allocate memory
  K <- length(rho)
  N <- nrow(theta)
  Kappa <- matrix(NA, nrow=N, ncol=K)

  # Compute inverses and log-determinants for each covariance matrix
  inv_sigma <- vector("list", length=K)
  for (k in 1:K){
    inv_sigma[[k]] <- solveV(Sigma[, , k], ldV=TRUE)
  }
  log_det <- sapply(inv_sigma, FUN=function(z){ z$ldV })
  log_rho <- log(rho)

  # Compute the Kappa matrix
  for (i in 1:N){
    for (k in 1:K){
      temp <- as.matrix(theta[i, ] - mu[k, ])
      Kappa[i, k] <- log_rho[k] - (t(temp) %*% inv_sigma[[k]]$y %*% temp + log_det[k]) / 2
    }
  }

  # Compute the matrix of multinomial probabilties
  Lambda <- exp(Kappa)
  Lambda <- Lambda / rowSums(Lambda)

  Z <- rcategorical(t(Lambda))

  return(list(z=Z, Lambda=Lambda))
}

## Helper functions for update functions ---------------------------------------

#' Sample from a categorical distribution.
#'
#'
#' Performs one draw from a categorical distribution (i.e., a multinomial distribution with size `n = 1`) for each of multiple probability vectors.
#'
#' @param prob An `n_cat x n_prob` matrix of probability vectors, each of which is a column of `prob`.  The entries of `prob` must be nonnegative, but are internally normalized such that each column sums to one.
#' @return A vector of length `n_prob` of integers between 1 and `n_cat` indicating the category draw for each column of `prob`.
#' @author Martin Lysy
rcategorical <- function(prob) {
  if(any(prob < 0)) stop("prob must contain only nonnegative entries.")
  cprob <- apply(prob, 2, cumsum)
  u <- runif(ncol(prob)) * cprob[nrow(prob),]

  apply(sweep(cprob, 2, u) >= 0, 2, which.max)
}

#' Random sampling from the Dirichlet distribution.
#'
#'
#' @param n Number of random draws.
#' @param alpha Weight parameter: a vector of nonnegative entries.
#' @return A matrix of size `n x length(alpha)` of which each row is a random draw.
#' @author Martin Lysy
rdirichlet <- function(n, alpha) {
  K <- length(alpha) # number of categories
  X <- matrix(rgamma(n*K, shape = alpha), K, n)

  drop(t(sweep(X, 2, colSums(X), "/")))
}

#' Solve method for variance matrices.
#'
#'
#' Calculates `y = V^{-1} x` where `V` is a variance matrix.
#'
#' @param V Variance matrix.
#' @param x Optional vector or matrix for which to solve system of equations.  If missing calculates inverse matrix.
#' @param ldV Optionally compute log determinant as well.
#'
#' @return Either a matrix or vector `y = V^{-1} x`, or if `ldV = TRUE`, a list with elements `y` and `ldV = log(det(V))`.
#' @details This function is faster and more stable than `base::solve()` when `V` is known to be positive-definite.
#' @author Martin Lysy
solveV <- function(V, x, ldV = FALSE) {
  C <- chol(V) # cholesky decomposition
  if(missing(x)) x <- diag(nrow(V))
  # solve is O(ncol(C) * ncol(x)) with triangular matrices
  # using backward subsitution
  ans <- backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
  if(ldV) {
    ldV <- 2 * sum(log(diag(C)))
    ans <- list(y = ans, ldV = ldV)
  }

  ans
}

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
#'
#' @references Arthur, D., Vassilvitskii, S. "k-means++: the advantages of careful seeding" *Proceedings of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms.* (2007): 1027–1035. <http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf>.
#' @author Martin Lysy
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
