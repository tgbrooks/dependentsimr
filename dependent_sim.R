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
  } else if (type == "DESeq2") {
    fit <- fit_deseq(data)
    marginals <- fit$marginals
    transformed_data <- fit$transformed_data
  } else {
    stop("'type' must be one of: 'normal', 'poisson', 'DESeq2'")
  }

  # Cov structure from PCA of the (transformed to have normal marginals) data
  udv <- svd(transformed_data, nu=rank, nv=0)
  udv$d <- udv$d / sqrt(dim(data)[2])

  # We compute variances without centering the data - since the transformation to normal
  # already 'centers' it (same reason we don't center/scale before SVD)
  variances <- apply(transformed_data^2, 1, mean)

  return(list(
    type = type,
    rank = rank,
    marginals = marginals,
    cov = udv,
    var = variances,
    n_features = dim(transformed_data)[1],
    transformed_data = transformed_data
  ))
}

draw_from_multivariate_corr <- function(random_structure, n_samples, size_factors=NULL) {
  # Samples from a multivariate distribution that has:
  # - each variable has the specified marginal distribution
  # - the same covariance in the top principal components observed in (normalized) data
  # 'size_factors' is used only when random_structure$type is DESeq2, in which case it is a vector
  # of length n_samples with the size factors (1 means average library size, 2 means twice the read depth, etc.)
  # By default size_factors is set to all 1s, so all libraries have approximately the same size
  
  k <- random_structure$rank
  marginals <- random_structure$marginals
  type <- random_structure$type
  n_features <- random_structure$n_features
  pc <- random_structure$cov
  
  # Draw from the multivariate normal distribution with the dependence structure of the pc
  # but done efficiently by transforming to a standard normal
  indep_draws <- matrix(rnorm(k*n_samples), c(k, n_samples))
  sdev <- diag(head(random_structure$cov$d, k), nrow=k)
  pc_draws <- pc$u %*% sdev %*% indep_draws
  
  # Add in the missing variance to match the actual data
  # by drawing independent data with the appropriate variance
  pc_var <- apply((pc$u %*% sdev)^2, 1, sum) # variance we already accounted for
  missing_var <- pmax(random_structure$var - pc_var, 0) # variance left to include as purely independent
  indep_draws <- matrix(rnorm(n_features*n_samples, sd=rep(sqrt(missing_var), n_samples)), c(n_features, n_samples))

  transformed_draws <- pc_draws + indep_draws

  # Adjust the draws to have the correct marginals
  if (type == "normal") {
    return(transformed_draws * sqrt(marginals$var) + marginals$mean)
  } else if (type == "poisson") {
    return(qpois(pnorm(transformed_draws), lambda = marginals$lambda))
  } else if (type == "DESeq2") {
    if (is.null(size_factors)) { size_factors <- rep(1, n_samples) }
    s <- matrix(size_factors, nrow=n_features, ncol=n_samples, byrow=TRUE)
    mu <- s * marginals$q
    draws <- qnbinom(
      pnorm(transformed_draws),
      mu = mu,
      size = 1/marginals$dispersion,
    )
    draws[is.na(draws)] <- 0 # NAs comes from the all-zero rows in the original data
    return(draws)
  } else {
    stop("'marginals$type' must be one of: 'normal', 'poisson', 'DESeq2'")
  }
}

remove_dependence <- function(random_structure) {
  # Returns another random structure that has 0 dependence
  # for comparison with the dependent version
  new_structure = random_structure
  new_structure$cov$d = rep(0, random_structure$rank)
  return(new_structure)
}

fit_deseq <- function(data) {
  # Use DESeq2 to fit the marginal distributions (negative binomial)
  dummy <- as.matrix(rep(1, dim(data)[2]))
  rownames(dummy) <- colnames(data)
  colnames(dummy) <- "dummy"

  dds <- DESeq2::DESeqDataSetFromMatrix(
    as.matrix(data),
    dummy,
    ~ 1
  )
  dds <- DESeq2::DESeq(dds)

  # Use the DESeq-fit model to transform the data to approximate normal
  # DESeq has the following model:
  # Each X_{ij} ~ NB(mu_ij, alpha_i)
  #      mu_ij = s_j q_{ij}
  #      log_2(q_{ij}) = x_j Beta_i
  # We use the trivial linear model so x_j = 1 and Beta_i is just an Intercept term
  # s_j is the sample-specific size factor to account for differences in read depths
  # alpha_i is the gene's dispersion parameter
  # perform the transformation via negative binomial -> probability -> normal distribution
  q <- 2^mcols(dds)$Intercept
  s <- dds$sizeFactor
  mu <- matrix(s, nrow=nrow(data), ncol=ncol(data), byrow=TRUE) * matrix(q, nrow=nrow(data), ncol=ncol(data))

  if (is.null(assays(dds)$replaceCounts)) {
    counts <- data
  } else {
    # DESeq2 has removed outliers, so we use that data
    counts <- assays(dds)$replaceCounts
  }
  transformed_data <- qnorm(pnbinom(
    counts,
    mu = mu,
    size = 1/dispersions(dds),
  ))

  # DESeq gives NAs where it can't fit. We replace those with zeros
  transformed_data[is.na(transformed_data)] <- 0
  return(list(
    marginals = list(
      q = q,
      dispersion = dispersions(dds)
    ),
    transformed_data = transformed_data
  ))
}
