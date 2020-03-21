require("mniw")

# Let %% denote matrix-vector or matrix-matrix multiplication

#' Samples a matrix, theta_i | A \ {theta_i}, from distribution N(G_i %% (y_i - mu_zi) - mu_zi,G_i %% V_i)
#' where G_i = Sigma_zi %% (V_i - Sigma_zi)^(-1), zi is a randomly sampled index term
#' 
#' @param y An 'n x 1' row vector denoting the features of a single observation.
#' @param V Observed fisher information
#' @param mu An 'n x 1' row vector denoting the prior means.
#' @param Sigma An 'n x n' matrix denoting the prior variances.
#' @return An 'n x 1' matrix that is sampled from mixture normal random effects model.
#' @details The vector that is returned is  `theta_i | A \ \{theta_i\}` which had been sampled from
#' ```
#' `N(G_i(y_i - mu_z_i) - mu_z_i,G_i V_i)``
#' where `G_i = Sigma_z_i(V_i - Sigma_z_i)^(-1)`, and  `z_i` is a randomly sampled index term
#' ```

OneStepSampleTheta <- function(y,V,mu,Sigma){
  theta_i <- rRxNorm(1,y,V,mu, Sigma)
  theta_i
}
  