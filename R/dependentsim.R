#' Compute structure of dependency from a given data
#'
#' @param datasets A list of data matrices that we will mimic. All datasets have samples in the columns and features (e.g., genes, proteins) in the rows. All datasets must have the same samples in corresponding columns.
#' @param rank Number of PCA components to approximate the dependence structure of.
#' @param types The marginal distribution types ('normal', 'poisson', 'DESeq2', or 'empirical'), as a list with entries corresponding to datasets. If just a single value is provided, then it is used for all datasets.
#'
#' @return A random structure element suitable for use with draw_from_multivariate_corr().
#' @export
#'
#' @importFrom stats runif
#' @importFrom stats qpois
#' @importFrom stats ppois
#' @importFrom stats rnorm
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats pnbinom
#' @importFrom stats qnbinom
#' @importFrom stats quantile
#' @importFrom stats ecdf
#'
#' @examples
#' # A small amount of RNA-seq read count data
#' d <- c(
#'   3563,      3839,      3407,      3745,      4063,      3089,      3219,      4577,      3769,
#'   0,         0,         0,         0,         0,         0,         0,         0,         0,
#'   0,         0,         0,         0,         0,         0,         0,         0,         0,
#'   3626,      2630,      5470,      4230,      4739,      4499,      5862,      6524,      3621,
#'   317,       310,       254,       449,       431,       463,       292,       325,       299,
#'   0,         0,         0,         0,         0,         0,         0,         0,         0,
#'   0,         0,         0,         0,         0,         0,         0,         0,         0,
#'   0,         0,         0,         0,         0,         0,         0,         0,         0,
#'   0,         0,         0,         0,         0,         0,         0,         0,         0,
#'   3,         2,         0,         1,         5,         3,         9,         1,         2,
#'   2,         1,         0,         2,         0,         0,         0,         0,         0,
#'   0,         0,         0,         0,         0,         0,         0,         0,         0,
#'   1,         0,         0,         0,         1,         0,         0,         0,         0,
#'   1,         1,         2,         0,         0,         1,         0,         0,         1,
#'   0,         0,         0,         0,         0,         0,         0,         0,         0,
#'   984,      1014,       777,      1211,      1174,      1096,       973,      1193,      1083,
#'   0,         0,         0,         0,         0,         0,         0,         0,         0,
#'   16,        15,        16,         9,        23,        10,        25,        10,         8,
#'   4473,      3954,      4476,      3965,      4606,      3788,      4084,      4316,      3658,
#'   0,         0,         0,         0,         0,         0,         0,         1,         0
#' )
#'
#' read_counts <- matrix(d, nrow=20, ncol=9, byrow=TRUE)
#'
#' # Simulate draws mimicking that data
#' rs_deseq <- get_random_structure(list(data=read_counts), rank=2, type="DESeq2")
#' draws_deseq <- draw_from_multivariate_corr(rs_deseq, n_samples=30)
get_random_structure <- function(datasets, rank, types="normal") {
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
      marginals[[dataname]] <- data
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

#' Draw random samples from the given random structure
#'
#' @param random_structure The structure of the random data to draw, computed by get_random_structure()
#' @param n_samples The number of samples to generate
#' @param size_factors Vector of length n_samples of size factors for DESeq2-type data. By default, all samples are given 1. Not used if no datasets has type DESeq2.
#'
#' @return A list of generated datasets. List names correspond to those of the datasets used to generate the random structure. Each dataset has n_samples columns.
#' @export
#'
#' @examples
#' #' # Generate example data with Sigma as its covariance matrix
#' library(MASS)
#' Sigma = matrix(c(
#'   1,      0.8,    0,  0,  0,  0,
#'   0.8,    1,      0,  0,  0,  0,
#'   0,      0,      1,  0,  0,  0,
#'   0,      0,      0,  1,  0,  0,
#'   0,      0,      0,  0,  1,  0.3,
#'   0,      0,      0,  0,  0.3, 1
#' ), nrow=6, ncol=6)
#' norm_data <- t(mvrnorm(n=20, mu=c(0,0,0,0,0,0), Sigma=Sigma))
#'
#' # Simulate draws mimicking that data
#' rs_normal <- get_random_structure(list(data=norm_data), rank=2, type="normal")
#' draws_normal <- draw_from_multivariate_corr(rs_normal, n_samples=30)
draw_from_multivariate_corr <- function(random_structure, n_samples, size_factors=NULL) {
  k <- random_structure$rank
  n_features <- random_structure$n_features
  pc <- random_structure$cov

  # Draw from the multivariate normal distribution with the dependence structure of the pc
  # but done efficiently by transforming to a standard normal
  indep_draws <- matrix(rnorm(k*n_samples), c(k, n_samples))
  sdev <- diag(random_structure$pc_factor_sizes[1:k], nrow=k)
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

#' Remove all dependence in a random structure
#'
#' @param random_structure A random structure from get_random_structure()
#'
#' @return The random structure with dependency removed, so all data generated from it will be independent.
#' @export
#'
#' @examples
#' library(MASS)
#' Sigma = matrix(c(
#'   1,      0.8,    0,  0,  0,  0,
#'   0.8,    1,      0,  0,  0,  0,
#'   0,      0,      1,  0,  0,  0,
#'   0,      0,      0,  1,  0,  0,
#'   0,      0,      0,  0,  1,  0.3,
#'   0,      0,      0,  0,  0.3, 1
#' ), nrow=6, ncol=6)
#' norm_data <- t(mvrnorm(n=20, mu=c(0,0,0,0,0,0), Sigma=Sigma))
#'
#' # Simulate draws mimicking that data but without any dependence
#' rs_normal <- get_random_structure(list(data=norm_data), rank=2, type="normal")
#' rs_indep <- remove_dependence(rs_normal)
#' draws_indep <- draw_from_multivariate_corr(rs_indep, n_samples=30)
remove_dependence <- function(random_structure) {
  new_structure = random_structure
  new_structure$pc_factor_sizes = rep(0, random_structure$rank)
  return(new_structure)
}

#' @importFrom stats pnbinom
#' @importFrom stats runif
#' @importFrom stats qnorm
fit_deseq <- function(data) {
  rlang::check_installed("DESeq2", reason="Using type='DESeq2' requires DESeq2")
  rlang::check_installed("S4Vectors", reason="Using type='DESeq2' requires S4Vectors")
  rlang::check_installed("SummarizedExperiment", reason="Using type='DESeq2' requires SummarizedExperiment")

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
  q <- 2^S4Vectors::mcols(dds)$Intercept
  s <- dds$sizeFactor
  mu <- matrix(s, nrow=nrow(data), ncol=ncol(data), byrow=TRUE) * matrix(q, nrow=nrow(data), ncol=ncol(data))

  if (is.null(SummarizedExperiment::assays(dds)$replaceCounts)) {
    counts <- data
  } else {
    # DESeq2 has removed outliers, so we use that data
    counts <- SummarizedExperiment::assays(dds)$replaceCounts
  }

  # this is a discrete distribution
  # So we 'smear' the probability of each bin out when converting to normal
  disp <- DESeq2::dispersions(dds)
  lower <- pnbinom(counts-1, mu = mu, size = 1/disp)
  upper <- pnbinom(counts, mu = mu, size = 1/disp)
  p <- matrix(runif(nrow(counts)*ncol(counts), lower, upper), nrow=nrow(counts), ncol=ncol(counts))

  transformed_data <- qnorm(p)

  # DESeq gives NAs where it can't fit. We replace those with zeros
  transformed_data[is.na(transformed_data)] <- 0
  return(list(
    marginals = list(
      q = q,
      dispersion = disp
    ),
    transformed_data = transformed_data
  ))
}
