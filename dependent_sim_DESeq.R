library(DESeq2)

get_random_structure_DEseq2 <- function(data, rank) {
  # Use DESeq2 to compute the marginal distributions of the data
  # Assuming all come from the same group (i.e., no predictor variables
  # in the DESeq2 model except for the constant 1).
  # Then determine a low-rank dependence structure appropriate for simulating
  # random data draws with draw_from_multivariate_corr_DESeq2
  # data is expected to be a matrix of unnormalized counts appropriate
  # for input to DESeq2.
  
  # Use DESeq2 to fit the marginal distributions (negative binomial)
  dummy <- as.matrix(rep(1, dim(data)[2]))
  rownames(dummy) <- colnames(data)
  colnames(dummy) <- "dummy"
  
  dds <- DESeqDataSetFromMatrix(
    as.matrix(data),
    dummy,
    ~ 1
  )
  dds <- DESeq(dds)
  
  # Use the DESeq-fit model to transform the data to approximate normal
  # DESeq has the following model:
  # Each X_{ij} ~ NB(mu_ij, alpha_i)
  #      mu_ij = s_j q_{ij}
  #      log_2(q_{ij}) = x_j Beta_i
  # We use the trivial linear model so x_j = 1 and Beta_i is just an Intercept term
  # s_j is the sample-specific size factor to account for differences in read depths
  # alpha_i is the gene's dispersion parameter
  # perform the transformation via negative binomial -> probability -> normal distribution
  transformed_data <- qnorm(pnbinom(
    data,
    mu = assays(dds)$mu,
    size = 1/dispersions(dds),
  ))
  
  # DESeq gives NAs where it can't fit. We replace those with zeros
  transformed_data[is.na(transformed_data)] <- 0

  # Cov structure from PCA of the (transformed to have normal marginals) data
  pc <- prcomp(t(transformed_data), rank. = rank, center = FALSE, scale. = FALSE)
  # We compute variances without centering the data - since the transformation to normal
  # already 'centers' it (same as center=FALSE in prcomp, otherwise we get negatives in the draw function)
  variances <- apply(transformed_data^2, 1, sum)

  return(list(
    type = "DESeq2",
    rank = rank,
    marginals = list(
      sizeFactor = dds$sizeFactor,
      q = assays(dds)$mu[,1] / dds$sizeFactor[1],
      dispersion = dispersions(dds)
    ),
    cov = pc,
    var = variances,
    n_features = dim(transformed_data)[1],
    transformed_data = transformed_data
  ))
}

draw_from_multivariate_corr_DESeq2 <- function(random_structure, n_samples, size_factors = NULL) {
  # Samples from a multivariate distribution that has:
  # - each variable has a negative binomial marginal distribution from the DESeq2 model
  # - the same covariance in the top k principal components observed in (normalized) data
  # - The size factors provided (constant 1s for all samples if omitted, giving identical library sizes for each)
  
  k <- random_structure$rank
  marginals <- random_structure$marginals
  type <- random_structure$type
  n_features <- random_structure$n_features
  pc <- random_structure$cov
  if (is.null(size_factors)) { size_factors <- rep(1, n_samples) }
  
  if (length(size_factors) != n_samples) { stop("size_factors must be of length n_samples")}
  if (random_structure$type != 'DESeq2') { stop("random structure was not generated via DESeq2")}
  
  # Draw from the multivariate normal distribution with the dependence structure of the pc
  # but done efficiently by transforming to a standard normal
  indep_draws <- matrix(rnorm(k*n_samples), c(k, n_samples))
  sdev <- diag(head(random_structure$cov$sdev, k), nrow=k)
  pc_draws <- pc$rotation %*% sdev %*% indep_draws
  
  # Add in the missing variance to match the actual data
  # by drawing independent data with the appropriate variance
  missing_var <- pmax(random_structure$var - (pc$rotation %*% sdev)^2, 0)
  #print(missing_var)
  indep_draws <- matrix(rnorm(n_features*n_samples, sd=rep(sqrt(missing_var), n_samples)), c(n_features, n_samples))
  
  # Random draws with normal marginals
  transformed_draws <- pc_draws + indep_draws
  
  # Correct the draws to have the correct marginals
  s <- matrix(size_factors, nrow=n_features, ncol=n_samples, byrow=TRUE)
  mu <- s * marginals$q
  draws <- qnbinom(
    pnorm(transformed_draws),
    mu = mu,
    size = 1/marginals$dispersion,
  )
  draws[is.na(draws)] <- 0 # NAs comes from the all-zero rows in the original data

  return(draws)
}
