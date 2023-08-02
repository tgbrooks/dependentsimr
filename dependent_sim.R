## NOTES
# - https://www.cse.unr.edu/~fredh/papers/conf/196-agpuantaafhdms/paper.pdf
#   GPU accelerated version of NORTA
# - https://citeseerx.ist.psu.edu/doc/10.1.1.48.281
#   Original(?) NORTA paper
# - DOI: 10.1002/asmb.901
#   For poisson variables specifically, gives efficient approximation to find the right input correlation matrix

X <- rbind( c(5, 4, 3), c(6, 10, 12), c(0, 0, 1), c(5, 4, 3.2), c(10, 8.2, 6.1))


covariance_for_simulation <- function(data, k) {
  # Finds a covariance matrix from the data that is appropriate for simulation
  # by assuming that all the independence lies in the top k principal components
  # but still makes the variance of each individual variable equal to the observed variance
  # The output is a covariance matrix, suitable for use with mvrnorm as the Sigma parameter
  
  # Compute the (sample) variances of each measurement in the data
  variances <- apply(data, c(1), var)
  
  # Perform principal component analysis
  pc <- prcomp(t(data), rank. = k)
  sim_cov <- pc$rotation %*% diag(head(pc$sdev, k), nrow=k)^2 %*%t(pc$rotation)
  diag(sim_cov) <- variances
  return(sim_cov)
}


sample_from_covariance_like <- function(data, k, n_samples) {
  # Samples from a multivariate normal distribution that has:
  # - each variable has the same variance as in the observed data
  # - the same variance in the top k principal components observed in data
  # The mean is assumed to be zero
  
  # Compute the (sample) variances of each measurement in the data
  variances <- apply(data, c(1), var)
  
  # Perform principal component analysis
  pc <- prcomp(t(data), rank. = k)
  
  # Draw from the multivariate normal distribution with the dependence structure of the pc
  # but done efficiently by transforming a standard normal
  indep_draws <- matrix(rnorm(k*n_samples), c(k, n_samples))
  sdev <- diag(head(pc$sdev, k), nrow=k)
  pc_draws <- pc$rotation %*% sdev %*% indep_draws
  
  # Add in the missing variance to match the actual data
  # by drawing independent data with the appropriate variance
  missing_var <- variances - (pc$rotation %*% sdev)^2
  n_features = dim(data)[1]
  indep_draws <- matrix(rnorm(n_features*n_samples, sd=rep(sqrt(missing_var), n_samples)), c(n_features, n_samples))

  return(pc_draws + indep_draws)
}

get_random_structure <- function(data, rank, type="normal") {
  # Compute the marginal and covariance structure from given data
  # Computes only the top 'rank' part of the covariance structure
  # type is one of 'normal', 'poisson' for the marginal distributions involved
  
  # First process the variables of the specified types
  # We compute each marginal distribution and transform the
  # data to have standard normal marginals
  if (type == "normal") {
    var <- apply(data, 1, var)
    mean <- apply(data, 1, mean)
    marginals <- list(
      var = var,
      mean = mean
    )
    # Trivial transform: just standardize the already normal data
    transformed_data <- (data - mean) / var
    
  } else if (type == "poisson") {
    marginals <- list(
      lambda = apply(data, 1, mean)
    )
    transformed_data <- qnorm(ppois(data, marginals$lambda))
  } else {
    stop("'type' must be one of: 'normal', 'poisson'")
  }
  
  # Cov structure from PCA of the (transformed to have normal marginals) data
  pc <- prcomp(t(transformed_data), rank. = rank)
  
  return(list(
    type = type,
    k = k,
    marginals = marginals,
    cov = pc,
    var = apply(transformed_data, 1, var),
    n_features = dim(transformed_data)[1],
    transformed_data = transformed_data
  ))
}

draw_from_multivariate_corr <- function(random_structure, n_samples) {
  # Samples from a multivariate distribution that has:
  # - each variable has the specified marginal distribution
  # - the same covariance in the top k principal components observed in data
  
  k <- random_structure$k
  marginals <- random_structure$marginals
  type <- random_structure$type
  n_features <- random_structure$n_features
  pc <- random_structure$cov
  
  # Draw from the multivariate normal distribution with the dependence structure of the pc
  # but done efficiently by transforming to a standard normal
  indep_draws <- matrix(rnorm(k*n_samples), c(k, n_samples))
  sdev <- diag(head(random_structure$cov$sdev, k), nrow=k)
  pc_draws <- pc$rotation %*% sdev %*% indep_draws
  
  # Add in the missing variance to match the actual data
  # by drawing independent data with the appropriate variance
  # Each variable has variance 1 since it has been transformed to standard normal
  missing_var <- random_structure$var - (pc$rotation %*% sdev)^2
  print(missing_var)
  indep_draws <- matrix(rnorm(n_features*n_samples, sd=rep(sqrt(missing_var), n_samples)), c(n_features, n_samples))

  # Random draws with normal marginals
  transformed_draws <- pc_draws + indep_draws

  # Correct the draws to have the correct marginals
  if (type == "normal") {
    return(transformed_draws * sqrt(marginals$var) + marginals$mean)
  } else if (type == "poisson") {
    return(qpois(pnorm(transformed_draws), lambda = marginals$lambda))
  } else {
    stop("'marginals$type' must be one of 'normal' or 'poisson'")
  }
}

 ## Test
library(MASS)
norm_data <- t(mvrnorm(n=20, mu=c(0,0, 0), Sigma=matrix(c(1, 0.8, 0, 0.8, 1, 0, 0, 0, 1), c(3,3))))
pois_data <- qpois(pnorm(norm_data), lambda = c(50, 1000, 2))

rs <- get_random_structure(pois_data, rank=1, type="poisson")
draws <- draw_from_multivariate_corr(rs, 30)
