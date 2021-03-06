% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update-functions.R
\name{update_mu}
\alias{update_mu}
\title{Conditional update for \code{mu_k} by randomly sampling from the Multivariate Normal model with mean theta and variance proportional to Sigma (see details).}
\usage{
update_mu(theta, old_mu, Sigma, z)
}
\arguments{
\item{theta}{Estimated random effects of the mixture-normal model. A \verb{n x p} matrix.}

\item{old_mu}{Current estimated mu. A \verb{K x p} matrix.}

\item{Sigma}{Variance matricies of theta. A \verb{p x p x K} array.}

\item{z}{Estimated classes of the y values. A \verb{n x 1} vector.}
}
\value{
A \verb{K x 1} vector of randomly sampled \code{mu_k} values from the distribution \code{mu_k ~ Normal(theta_k, Sigma_k/N_k)}\preformatted{where `theta_k` is the mean of all `theta_i`s where `z_i = k` and `N_k` is the number of `z_i = k` (sum of indicator `z_k = k`).
}
}
\description{
Conditional update for \code{mu_k} by randomly sampling from the Multivariate Normal model with mean theta and variance proportional to Sigma (see details).
}
\details{
\code{K} is the number of classes. \code{p} is the dimension of one observations of \code{theta_i}, i.e. the number of features. \code{n} is the number of observations. If a class \code{k} does not appear in the z vector, the \code{k}th row of the return matrix is the \code{k}th row of currMu. If a class has no observations in the z vector this mu is not updated because having a count of zero causes a division by zero.
}
