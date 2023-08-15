library(DESeq2)

get_random_structure <- function(datasets, rank, types="normal") {
  # Compute the marginal and covariance structure from given data
  # Computes only the top 'rank' part of the covariance structure
  
  # datasets is a list of data matrices, with each entry a matrix, which all have common columns
  # (i.e., they all have n columns and the ith column of one corresponds to the ith column of another).
  # types is a list with values that are all in: 'normal', 'poisson', 'DESeq2', or 'empirical' specifying the
  # maginal distribution type to use for the corresponding data in datasets.
  # If just a single value is provided, then this is used for all data in datasets.

  if (length(types) == 1) {
    # use this type for all data matrices
    types <- lapply(datasets, function(x) types[1])
  }

  # datasets and types must have corresponding entries
  stopifnot(all(names(datasets) %in% names(types)))

  # all datasets must all have the same number of samples (columns)
  numcols <- lapply(datasets, ncol)
  stopifnot(length(unique(numcols)) == 1)
  n <- unique(numcols)[[1]]

  marginals <- list()
  transformed_data <- list()
  for (dataname in names(datasets)) {
    data <- datasets[[dataname]]
    type <- types[[dataname]]

    # First process the variables of the specified types
    # We compute each marginal distribution and transform the
    # data to have standard normal marginals
    if (type == "normal") {
      var <- apply(data, 1, var)
      mean <- apply(data, 1, mean)
      marginals[[dataname]] <- list(
        var = var,
        mean = mean
      )
      # Trivial transform: just standardize the already normal data
      transformed_data[[dataname]] <- (data - mean) / var
    } else if (type == "poisson") {
      lambda <- apply(data, 1, mean)
      marginals[[dataname]] <- list(
        lambda = lambda
      )
      transformed_data[[dataname]] <- qnorm(ppois(data, lambda))
    } else if (type == "DESeq2") {
      fit <- fit_deseq(data)
      marginals[[dataname]] <- fit$marginals
      transformed_data[[dataname]] <- fit$transformed_data
    } else if (type == "empirical") {
      marginals[dataname] <- data
      # Transform by the empirical CDF, with mass 'smeared' out so that the full range of 0-1 is possible
      lower <- t(apply(data, 1, function(x) ecdf(x)(x-1e-8)))
      upper <- t(apply(data, 1, function(x) ecdf(x)(x)))
      p <- matrix(runif(nrow(data)*ncol(data), lower, upper), nrow=nrow(data), ncol=ncol(data))
      transformed_data[[dataname]] <- qnorm(p) # normalize
    } else {
      stop("'type' must be one of: 'normal', 'poisson', 'DESeq2', 'empirical'")
    }
  }

  transformed_data_matrix <- do.call(rbind, transformed_data)

  # Cov structure from PCA of the (transformed to have normal marginals) data
  udv <- svd(transformed_data_matrix, nu=rank, nv=0)
  udv$d <- udv$d / sqrt(n-1)

  # We compute variances without centering the data - since the transformation to normal
  # already 'centers' it (same reason we don't center/scale before SVD)
  # Require all to have at least variance 1 - we expect them all to have exactly variance 1
  # since we have transformed to standard normal, but they do not exactly equal 1 in general
  variances <- pmax(apply(transformed_data_matrix^2, 1, sum) / (n-1), 1)

  # Compute the appropriate scaling for each the PCs, accounting for the overlap
  # with the independent data that gets added on top
  # In real omics data, this is very close to just doing: pc_factor_sizes <- udv$d
  # but for smaller datasets this has a bigger impact
  rhs <- vapply(1:rank, function(k) {
    udv$d[k]^2 - sum(udv$u[ ,k]^2 * variances)
  }, FUN.VALUE=1.0)
  U2 <- udv$u^2
  lhs <- diag(rank) - t(U2) %*% U2
  pc_factor_sizes <- sqrt(solve(lhs, rhs))

  return(list(
    type = types,
    rank = rank,
    marginals = marginals,
    cov = udv,
    pc_factor_sizes = pc_factor_sizes,
    var = variances,
    n_features = dim(transformed_data_matrix)[1],
    n_features_list = lapply(datasets, nrow),
    n_samples = n,
    transformed_data = transformed_data,
    rownames = lapply(datasets, rownames)
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
  n_features <- random_structure$n_features
  pc <- random_structure$cov
  
  # Draw from the multivariate normal distribution with the dependence structure of the pc
  # but done efficiently by transforming to a standard normal
  indep_draws <- matrix(rnorm(k*n_samples), c(k, n_samples))
  sdev <- diag(head(random_structure$pc_factor_sizes, k), nrow=k)
  pc_draws <- pc$u %*% sdev %*% indep_draws
  
  # Add in the missing variance to match the actual data
  # by drawing independent data with the appropriate variance
  pc_var <- apply((pc$u %*% sdev)^2, 1, sum) # variance we already accounted for
  missing_var <- pmax(random_structure$var - pc_var, 0) # variance left to include as purely independent
  indep_draws <- matrix(rnorm(n_features*n_samples, sd=rep(sqrt(missing_var), n_samples)), c(n_features, n_samples))

  transformed_draws_all <- pc_draws + indep_draws

  start_row <- 1
  draws <- list()
  for (dataname in names(random_structure$n_features_list)) {
    # Extract each dataset out
    n_features <- random_structure$n_features_list[[dataname]]
    end_row <- start_row + n_features - 1
    transformed_draws <- transformed_draws_all[start_row:end_row,]
    marginals <- random_structure$marginals[[dataname]]
    type <- random_structure$type[[dataname]]

    # Adjust the normal draws to have the correct marginals
    if (type == "normal") {
      draws[[dataname]] <- transformed_draws * sqrt(marginals$var) + marginals$mean
    } else if (type == "poisson") {
      draws[[dataname]] <- qpois(pnorm(transformed_draws), lambda = marginals$lambda)
    } else if (type == "DESeq2") {
      #TODO: can we support multiple size-factors per sample in case we have multiple DESeq2 datasets simultaneously?
      if (is.null(size_factors)) { size_factors <- rep(1, n_samples) }
      s <- matrix(size_factors, nrow=n_features, ncol=n_samples, byrow=TRUE)
      mu <- s * marginals$q
      d <- qnbinom(
        pnorm(transformed_draws),
        mu = mu,
        size = 1/marginals$dispersion,
      )
      d[is.na(d)] <- 0 # NAs comes from the all-zero rows in the original data
      draws[[dataname]] <- d
    } else if (type == "empirical") {
      p <- pnorm(transformed_draws)
      d <- matrix(NA, nrow=nrow(transformed_draws), ncol=ncol(transformed_draws))
      for (i in 1:nrow(d)) {
        d[i,] <- quantile(marginals[i,], p[i,], type=1)
      }
      draws[[dataname]] <- d
    } else {
      stop("'marginals$type' must be one of: 'normal', 'poisson', 'DESeq2', 'empirical'")
    }
    rownames(draws[[dataname]]) <- random_structure$rownames[[dataname]]
    
    start_row <- start_row + n_features
  }

  return(draws)
}

remove_dependence <- function(random_structure) {
  # Returns another random structure that has 0 dependence
  # for comparison with the dependent version
  new_structure = random_structure
  new_structure$pc_factor_sizes = rep(0, random_structure$rank)
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

  # this is a discrete distribution
  # So we 'smear' the probability of each bin out when converting to normal
  lower <- pnbinom(counts-1, mu = mu, size = 1/dispersions(dds))
  upper <- pnbinom(counts, mu = mu, size = 1/dispersions(dds))
  p <- matrix(runif(nrow(counts)*ncol(counts), lower, upper), nrow=nrow(counts), ncol=ncol(counts))

  transformed_data <- qnorm(p)

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
