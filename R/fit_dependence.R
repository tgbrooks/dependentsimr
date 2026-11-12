DEPENDENCE_METHODS <- c("pca", "spiked Wishart", "corpcor", "independent")

#' Fit the dependence between features
#'
#' The second step of a simulation: estimate the dependence between features
#' (genes, proteins, etc.) from the residuals of fit marginals, as the covariance
#' matrix of a Gaussian copula. Since the residuals account for the design given
#' to [fit_marginals()], the dependence is that remaining after removing the
#' effects of known covariates. The covariance is represented as low-rank plus
#' diagonal so that it scales to omics-sized data.
#'
#' Methods:
#' * `"pca"`: dependence from the top `rank` principal components of the residuals.
#'   Small values of `rank` (e.g., 2 or 3) are typical unless there are many samples.
#' * `"spiked Wishart"`: also uses the top `rank` components, but fits their sizes
#'   with a spiked Wishart model to correct for the bias in sample eigenvalues.
#'   `rank` can be up to the number of samples minus the number of design columns.
#'   Requires the `sparsesvd` and `Matrix` packages.
#' * `"corpcor"`: the shrinkage covariance estimator of the `corpcor` package. Takes
#'   no parameters but often underestimates dependence. Requires `corpcor`.
#' * `"independent"`: no dependence between features, for comparison.
#'
#' @param marginals Fit marginals from [fit_marginals()]
#' @param method One of `"pca"`, `"spiked Wishart"`, `"corpcor"`, or `"independent"`
#' @param rank Number of components of dependence. Used by `"pca"` and `"spiked Wishart"`.
#'
#' @return An object of class `dependentsimr_dependence` for use with [simulate_data()].
#' @export
#'
#' @examples
#' true_means <- c(0, 1, 10, 1000)
#' data <- rpois(4*12, true_means) |> matrix(4, 12)
#' marginals <- fit_marginals(data, method = "poisson")
#'
#' # Try out different dependence methods without refitting the marginals
#' dep_pca <- fit_dependence(marginals, method = "pca", rank = 2)
#' dep_indep <- fit_dependence(marginals, method = "independent")
fit_dependence <- function(marginals, method, rank = 2) {
  check_marginals(marginals)
  if (!(length(method) == 1 && method %in% DEPENDENCE_METHODS)) {
    stop("Dependence 'method' must be one of: ", paste0("'", DEPENDENCE_METHODS, "'", collapse = ", "))
  }

  residuals <- do.call(rbind, lapply(marginals$datasets, function(ds) ds$residuals))
  constant <- unlist(lapply(marginals$datasets, function(ds) ds$constant), use.names = FALSE)
  m <- nrow(residuals)

  # The residuals are constrained to the space orthogonal to the design, so they are
  # not independent across samples. Rotating onto an orthonormal basis of that space
  # gives n - p columns that are (approximately) iid draws of the features' dependence.
  # For an intercept-only design this is equivalent to centering and dividing by n-1.
  design <- marginals$design
  n <- nrow(design)
  p <- ncol(design)
  basis <- qr.Q(qr(design), complete = TRUE)[, (p + 1):n, drop = FALSE]
  Z <- residuals %*% basis
  dimnames(Z) <- NULL
  n_eff <- n - p

  # We compute variances without centering - the transformation to normal already
  # 'centers' the data. Require all to have at least variance 1 - we expect them all
  # to have variance 1 since we have transformed to standard normal, but they do not
  # exactly equal 1 in general
  variances <- pmax(rowSums(Z^2) / n_eff, 1)

  result <- list(
    method = method,
    n_features = m,
    n_features_list = lapply(marginals$datasets, function(ds) ds$n_features),
    rownames = lapply(marginals$datasets, function(ds) ds$rownames),
    n_eff = n_eff
  )

  if (method %in% c("pca", "spiked Wishart")) {
    if (!(is.numeric(rank) && length(rank) == 1 && rank >= 1 && rank == round(rank))) {
      stop("'rank' must be a positive integer")
    }
    if (rank > min(n_eff, m)) {
      stop(
        "'rank' can be at most the number of samples minus the number of design columns (",
        n_eff, ") and at most the number of features (", m, ")"
      )
    }
    result$rank <- rank
  }

  if (method == "pca") {
    # Cov structure from PCA of the (transformed to have normal marginals) data
    udv <- svd(Z, nu = rank, nv = 0)
    udv$d <- udv$d / sqrt(n_eff)

    # Compute the appropriate scaling for each the PCs, accounting for the overlap
    # with the independent data that gets added on top
    # In real omics data, this is very close to just doing: pc_factor_sizes <- udv$d
    # but for smaller datasets this has a bigger impact
    rhs <- vapply(1:rank, function(k) {
      udv$d[k]^2 - sum(udv$u[, k]^2 * variances)
    }, FUN.VALUE = 1.0)
    U2 <- udv$u^2
    lhs <- diag(rank) - t(U2) %*% U2
    result$cov <- udv
    result$pc_factor_sizes <- sqrt(solve(lhs, rhs))
    result$var <- variances
  } else if (method == "spiked Wishart") {
    rlang::check_installed("sparsesvd", reason = "Using method='spiked Wishart' requires sparsesvd package")
    rlang::check_installed("Matrix", reason = "Using method='spiked Wishart' requires Matrix package")

    udv <- svd(Z, nu = rank, nv = 0)

    # Don't count zero-variance variables as variables for this
    num_vars <- m - sum(constant)
    num_eigs <- min(n_eff, num_vars)
    fit <- match_with_spiked_wishart(
      desired_eigenvalues = udv$d[1:num_eigs],
      rank = rank,
      num_observations = n_eff,
      num_variables = num_vars,
      num_iterations = 20,
      num_samples_per_iter = 300
    )
    result$cov <- udv
    result$spiked_sd <- fit$spiked_sd
    result$population_sd <- fit$population_sd
    result$var <- variances
  } else if (method == "corpcor") {
    rlang::check_installed("corpcor", reason = "Using method='corpcor' requires corpcor package")
    if (n_eff < 2) {
      stop("method = 'corpcor' requires the number of samples to exceed the number of design columns by at least 2")
    }

    # Emulate the estimator from corpcor's cov.shrink
    # using corpcor's estimates for the two lambdas
    # This shrinks the correlation matrix to the identity and also shrinks the variances towards the median variance.
    # The end result matrix can be expressed as a decomposition into a independent (i.e., diagonal)
    # part and dependent (non-diagonal) but low-rank part. In particular, we just calculate the square root of the
    # dependent part which is only m x n_eff instead of m x m
    # NOTE: we've already normalized each gene, so the variance shrinkage doesn't do much but we include it anyway
    # NOTE: we zero out the variables that zero variance since otherwise they force the shrinkage of all correlation towards 0
    dat <- Z
    dat[constant, ] <- 0
    v <- rowSums(dat^2) / n_eff # variances (uncentered, see above)

    # corpcor centers the data internally, which would remove another degree of freedom.
    # So we give it the data embedded into n_eff+1 columns that are already centered:
    # this has the same cross-products as 'dat' so centering is a no-op.
    centering <- t(qr.Q(qr(matrix(1, n_eff + 1, 1)), complete = TRUE)[, -1, drop = FALSE])
    # All-zero rows contribute nothing to corpcor's lambda estimates but do trigger warnings
    informative <- rowSums(dat^2) > 0
    dat_centered <- dat[informative, , drop = FALSE] %*% centering
    lambda <- corpcor::estimate.lambda(t(dat_centered), verbose = FALSE)
    lambda.var <- corpcor::estimate.lambda.var(t(dat_centered), verbose = FALSE)

    med_v <- stats::median(v) # we shrink towards the median variance
    indep_part <- lambda * (lambda.var * med_v + (1 - lambda.var) * v)
    dep_part <- sqrt((1 - lambda) / n_eff) * sweep(dat, 1, sqrt(lambda.var * med_v / v + (1 - lambda.var)), "*")
    dep_part[is.nan(dep_part)] <- 0 # NaNs come in the 0 variance features
    # The shrunk covariance is then: diag(indep_part) + dep_part %*% t(dep_part)
    result$indep_part <- indep_part
    result$dep_part <- dep_part
    result$var <- variances
  }

  structure(result, class = "dependentsimr_dependence")
}

#' @export
print.dependentsimr_dependence <- function(x, ...) {
  cat("dependentsimr dependence: method '", x$method, "'", sep = "")
  if (!is.null(x$rank)) cat(", rank", x$rank)
  cat(",", x$n_features, "features\n")
  invisible(x)
}

check_dependence <- function(dependence) {
  if (!inherits(dependence, "dependentsimr_dependence")) {
    stop("'dependence' must be the result of fit_dependence()")
  }
}

# Draw from the Gaussian copula: a n_features x n_samples matrix with
# approximately standard normal marginals and the fit dependence
draw_latent <- function(dependence, n_samples) {
  m <- dependence$n_features
  method <- dependence$method
  if (method == "pca" || method == "spiked Wishart") {
    k <- dependence$rank
    pc <- dependence$cov
    sizes <- if (method == "pca") dependence$pc_factor_sizes else dependence$spiked_sd

    # Draw from the multivariate normal distribution with the dependence structure of the pc
    # but done efficiently by transforming to a standard normal
    sdev <- diag(sizes, nrow = k)
    pc_draws <- pc$u %*% sdev %*% matrix(rnorm(k * n_samples), k, n_samples)

    # Add in the missing variance to match the actual data
    # by drawing independent data with the appropriate variance
    if (method == "pca") {
      pc_var <- rowSums((pc$u %*% sdev)^2) # variance we already accounted for
      missing_var <- pmax(dependence$var - pc_var, 0) # variance left to include as purely independent
    } else {
      missing_var <- abs(dependence$population_sd)
    }
    indep_draws <- matrix(rnorm(m * n_samples, sd = rep(sqrt(missing_var), n_samples)), m, n_samples)
    pc_draws + indep_draws
  } else if (method == "corpcor") {
    k <- ncol(dependence$dep_part)
    dep_draws <- dependence$dep_part %*% matrix(rnorm(k * n_samples), k, n_samples)
    indep_draws <- matrix(rnorm(m * n_samples, sd = rep(sqrt(dependence$indep_part), n_samples)), m, n_samples)
    dep_draws + indep_draws
  } else if (method == "independent") {
    matrix(rnorm(m * n_samples), m, n_samples)
  } else {
    stop("Unknown dependence method '", method, "'")
  }
}
