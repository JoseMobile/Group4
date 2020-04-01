require("mniw")

# Let %% denote matrix-vector or matrix-matrix multiplication

#' Samples a matrix, theta_i | A \ {theta_i}, from distribution N(G_i %% (y_i - mu_zi) - mu_zi,G_i %% V_i)
#' where G_i = Sigma_zi %% (V_i - Sigma_zi)^(-1), zi is a randomly sampled index term
#' 
#' @param Y An matrix denoting the features of a single observation.
#' @param V An array of Observed fisher information matrices
#' @param mu A matrix denoting the prior column means of each cluster.
#' @param Sigma An array of 'n x n' matrices denoting the prior variances for each cluster.
#' @return An array 'n x 1' vectors that is sampled from mixture normal random effects model.
#' @details The vector that is returned is  `theta_i | A \ \{theta_i\}` which had been sampled from
#' ```
#' `N(G_i(y_i - mu_z_i) - mu_z_i,G_i V_i)``
#' where `G_i = Sigma_z_i(V_i - Sigma_z_i)^(-1)`, and  `z_i` is a randomly sampled index term
#' ```

OneStepSampleTheta <- function(Y,V,mu,Sigma){
  theta <- rRxNorm(1,Y,V,mu, Sigma)
  theta
}
  