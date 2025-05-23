% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spiked_wishart.R
\name{sample_spiked_wishart_and_jac}
\alias{sample_spiked_wishart_and_jac}
\title{Efficiently sample the eigenvalues of a random Wishart matrix with spiked eigenvalues and compute the derivative
of the singular values with respect to spiked eigenvalues}
\usage{
sample_spiked_wishart_and_jac(
  spiked_sd,
  num_observations,
  num_variables,
  population_sd = 1,
  num_eigs = 0
)
}
\arguments{
\item{spiked_sd}{The spiked standard deviations}

\item{num_observations}{The number of observations (aka samples or columns)}

\item{num_variables}{The number of variables (aka features or rows)}

\item{population_sd}{the standard deviation of all non-spiked components (num_variables - length(spiked_sd) of them)}

\item{num_eigs}{The number of eigenvalues to compute. If 0 compute all of them using dense matrix routines. If greater than zero, use sparse matrices and compute that many top eigenvalues.}
}
\value{
List with both a vector of random singular values of G where G is a random num_variables x num_observations matrix with iid columns from N(0, Sigma) where Sigma is diagonal with entries spiked_sd and all the remaining are population_sd.
and also the Jacobian, where \verb{[[i,j]]} is the derivative of the ith singular value with respect to the jth spiked SD
}
\description{
Efficiently sample the eigenvalues of a random Wishart matrix with spiked eigenvalues and compute the derivative
of the singular values with respect to spiked eigenvalues
}
\examples{
set.seed(0)
# Sample eigenvalues of a covariance matrix of a 10 sample study with 1000 variables such that
# the top two (underlying true distribution of the data, not sample) principal components
# have SDs of 100 and 10 and the remaining 98 have 1
res = sample_spiked_wishart_and_jac(
  spiked_sd = c(500, 100),
  num_variables = 1000,
  num_observations = 10-1,
  num_eigs = 3
)
res$singular_vals # singular values (of G, i.e., square roots of the eigenvalues of W = G G^T)
res$jacobian # jacobian of the singular values (sqrt of the eigenvalues)
             # with respect to each of the spiked_sd's
res$pop_sd_grad # Gradient of the singular values with respect to the population_sd parameter
}
