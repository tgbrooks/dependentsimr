MARGINAL_METHODS <- c("normal", "poisson", "DESeq2", "empirical")

#' Fit the marginal distribution of each feature
#'
#' The first step of a simulation: fit a model for each feature (gene, protein,
#' etc.) individually, and compute residuals of the data under that model.
#' The residuals are randomized quantile residuals, which are standard normal
#' if the model is correct, and are what [fit_dependence()] uses to estimate the
#' dependence between features.
#'
#' Each feature's distribution may depend upon known sample-level covariates
#' through a design matrix. The models are:
#' * `"normal"`: \eqn{Y_{ij} \sim N(x_j \beta_i, \sigma_i^2)}, fit with [stats::lm.fit()].
#'   Coefficients are on the scale of the data.
#' * `"poisson"`: \eqn{Y_{ij} \sim Pois(\exp(x_j \beta_i))}, fit with [stats::glm.fit()].
#'   Coefficients are on the natural-log scale.
#' * `"DESeq2"`: \eqn{Y_{ij} \sim NB(s_j 2^{x_j \beta_i}, \alpha_i)} with size factors \eqn{s_j},
#'   fit with `DESeq2::DESeq()`. Coefficients are log2 fold changes, as in DESeq2.
#'   The dispersion \eqn{\alpha_i} uses DESeq2's gene-wise (unshrunk) estimate,
#'   since the shrunken estimates give unrealistic simulations.
#' * `"empirical"`: each feature takes the values observed in the data. Does not support a design.
#'
#' Since fitting (particularly with DESeq2) is the slowest step, the result can
#' be reused with multiple calls to [fit_dependence()].
#'
#' @param data The data to mimic, with features in rows and samples in columns.
#'   Either a single matrix or data frame, or a named list of them for multiple
#'   datasets (e.g., transcriptomics and proteomics) measured on the same samples
#'   (in corresponding columns). A data frame must be entirely numeric except
#'   optionally its first column, which is then used as the feature names.
#' @param method The marginal distribution: one of `"normal"`, `"poisson"`,
#'   `"DESeq2"`, or `"empirical"`. For multiple datasets, either a single value
#'   used for all of them or a named list with an entry for each dataset.
#' @param design Optional design matrix of known sample covariates, with one row
#'   per sample (in the same order as the columns of `data`) and one column per
#'   covariate, usually from [stats::model.matrix()]. If `NULL` (the default),
#'   samples are treated as identically distributed (an intercept-only model).
#'   If provided, [simulate_data()] must also be given a `new_design`.
#'
#' @return An object of class `dependentsimr_marginals`. Its `datasets` entry is a
#'   list with an element per dataset, each containing the `method`, the fit
#'   `coefficients` (features x covariates matrix, `NULL` for empirical),
#'   per-feature parameters (`sigma2` for normal, `dispersion` for DESeq2,
#'   `values` for empirical), and the `residuals` (features x samples).
#'   Use [add_effect()] to modify the coefficients for simulating differential expression.
#' @export
#'
#' @examples
#' # Simple data to mimic
#' true_means <- c(0, 1, 10, 1000)
#' data <- rpois(4*12, true_means) |> matrix(4, 12)
#'
#' marginals <- fit_marginals(data, method = "poisson")
#' dependence <- fit_dependence(marginals, method = "pca", rank = 2)
#' sim <- simulate_data(marginals, dependence, n_samples = 30)
#'
#' # With a known two-group structure in the data
#' group <- factor(rep(c("A", "B"), each = 6))
#' design <- model.matrix(~ group)
#' marginals <- fit_marginals(data, method = "poisson", design = design)
fit_marginals <- function(data, method, design = NULL) {
  input <- as_dataset_list(data)
  datasets <- input$datasets

  if (is.list(method)) {
    if (is.null(names(method)) || !all(names(datasets) %in% names(method))) {
      stop("When 'method' is a list, it must be named with an entry for each dataset")
    }
    methods <- method[names(datasets)]
  } else {
    if (length(method) != 1) {
      stop("'method' must be a single value or a named list with an entry per dataset")
    }
    methods <- lapply(datasets, function(x) method)
  }
  for (m in methods) {
    if (!(length(m) == 1 && m %in% MARGINAL_METHODS)) {
      stop("Marginal 'method' must be one of: ", paste0("'", MARGINAL_METHODS, "'", collapse = ", "))
    }
  }

  n <- unique(vapply(datasets, ncol, 1L))
  if (length(n) != 1) {
    stop("All datasets must have the same number of samples (columns)")
  }

  design_supplied <- !is.null(design)
  design <- check_design(design, n, sample_names = colnames(datasets[[1]]))
  if (design_supplied && any(methods == "empirical")) {
    stop("method = 'empirical' does not support a design matrix")
  }
  hat <- rowSums(qr.Q(qr(design))^2)
  if (any(hat > 1 - 1e-6)) {
    warning(
      "Some samples are perfectly fit by the design (e.g., they are alone in their group) ",
      "so their residuals are zero and they carry no information about dependence"
    )
  }

  fits <- list()
  for (dataname in names(datasets)) {
    Y <- datasets[[dataname]]
    fit <- switch(methods[[dataname]],
      normal = fit_normal(Y, design),
      poisson = fit_poisson(Y, design),
      DESeq2 = fit_deseq(Y, design),
      empirical = fit_empirical(Y)
    )
    fit$method <- methods[[dataname]]
    fit$constant <- apply(Y, 1, var) == 0
    fit$rownames <- rownames(Y)
    fit$n_features <- nrow(Y)
    rownames(fit$residuals) <- rownames(Y)
    colnames(fit$residuals) <- colnames(Y)
    if (!is.null(fit$coefficients)) {
      rownames(fit$coefficients) <- rownames(Y)
      colnames(fit$coefficients) <- colnames(design)
    }
    # Order fields consistently for printing/inspection
    first <- c("method", "coefficients")
    fits[[dataname]] <- fit[c(intersect(first, names(fit)), setdiff(names(fit), first))]
  }

  structure(
    list(
      datasets = fits,
      design = design,
      design_supplied = design_supplied,
      new_design_required = design_supplied,
      single = input$single,
      n_samples = n
    ),
    class = "dependentsimr_marginals"
  )
}

#' Add or set a covariate effect in fit marginals
#'
#' Adds a new coefficient (column of the design) to fit marginals, typically to
#' simulate differential expression between groups that were not present in the
#' original data. The simulated data then depends on the corresponding column of
#' `new_design` given to [simulate_data()], which becomes required.
#'
#' Note that the coefficients estimated from real data generally should not be
#' used as-is as the 'true' effects for a simulation: estimation noise makes every
#' feature appear to have some effect, inflating the amount of differential
#' expression compared to the original data. Instead, choose effects directly,
#' e.g. a random subset of genes with random log fold changes, and use
#' `overwrite = TRUE` if replacing a fit coefficient.
#'
#' @param marginals Fit marginals from [fit_marginals()]
#' @param name Name of the coefficient. Must match a column name of `new_design` in [simulate_data()].
#' @param effects The effect sizes, on the coefficient scale of the dataset's method:
#'   log2 fold change for `"DESeq2"`, natural log fold change for `"poisson"`, or
#'   additive shift for `"normal"`. Either a numeric vector with one entry per feature,
#'   or a named vector of effects for some features (matched by feature name, all others get 0).
#' @param dataset Name of the dataset to apply the effects to. May be omitted if there is only one.
#'   Other datasets get an effect of zero.
#' @param overwrite If `TRUE`, allow replacing an existing coefficient (such as one from the design).
#'
#' @return The modified marginals.
#' @export
#'
#' @examples
#' data <- rpois(100*12, rep(c(5, 50, 500, 5000), 25)) |> matrix(100, 12)
#' rownames(data) <- paste0("gene", 1:100)
#' marginals <- fit_marginals(data, method = "poisson")
#' dependence <- fit_dependence(marginals, method = "pca", rank = 2)
#'
#' # Ten genes get a log fold change of +/- 1
#' de_genes <- sample(rownames(data), 10)
#' lfc <- setNames(sample(c(-1, 1), 10, replace = TRUE), de_genes)
#' marginals_de <- add_effect(marginals, "treatment", lfc)
#'
#' new_design <- cbind("(Intercept)" = 1, treatment = rep(c(0, 1), each = 5))
#' sim <- simulate_data(marginals_de, dependence, new_design = new_design)
add_effect <- function(marginals, name, effects, dataset = NULL, overwrite = FALSE) {
  check_marginals(marginals)
  if (!(is.character(name) && length(name) == 1 && !is.na(name) && nzchar(name))) {
    stop("'name' must be a single string")
  }
  if (!is.numeric(effects)) {
    stop("'effects' must be numeric")
  }
  datanames <- names(marginals$datasets)
  if (is.null(dataset)) {
    if (length(datanames) > 1) {
      stop("'dataset' must be given when there are multiple datasets")
    }
    dataset <- datanames[[1]]
  }
  if (!(dataset %in% datanames)) {
    stop("No dataset named '", dataset, "'")
  }
  target <- marginals$datasets[[dataset]]
  if (target$method == "empirical") {
    stop("Cannot add an effect to an 'empirical' dataset")
  }

  n_features <- target$n_features
  if (!is.null(names(effects))) {
    if (is.null(target$rownames)) {
      stop("'effects' is named but the dataset has no feature names; give an unnamed vector with an entry per feature")
    }
    unknown <- setdiff(names(effects), target$rownames)
    if (length(unknown) > 0) {
      stop("Unknown feature names in 'effects': ", paste(utils::head(unknown, 5), collapse = ", "))
    }
    full <- rep(0, n_features)
    full[match(names(effects), target$rownames)] <- effects
    effects <- full
  } else if (length(effects) != n_features) {
    stop("'effects' must have one entry per feature (", n_features, ") or be named by feature")
  }

  existing <- coefficient_names(marginals)
  if (name %in% existing && !overwrite) {
    stop("Coefficient '", name, "' already exists; use overwrite = TRUE to replace it")
  }

  for (dataname in datanames) {
    ds <- marginals$datasets[[dataname]]
    if (is.null(ds$coefficients)) next # empirical: does not depend on the design
    value <- if (dataname == dataset) effects else 0
    if (name %in% colnames(ds$coefficients)) {
      ds$coefficients[, name] <- value
    } else {
      ds$coefficients <- cbind(ds$coefficients, value)
      colnames(ds$coefficients)[ncol(ds$coefficients)] <- name
    }
    marginals$datasets[[dataname]] <- ds
  }
  marginals$new_design_required <- TRUE
  marginals
}

#' @export
print.dependentsimr_marginals <- function(x, ...) {
  cat("dependentsimr marginals:", x$n_samples, "samples\n")
  for (dataname in names(x$datasets)) {
    ds <- x$datasets[[dataname]]
    cat("  ", dataname, ": ", ds$n_features, " features, method '", ds$method, "'\n", sep = "")
  }
  if (x$design_supplied) {
    cat("  design:", paste(colnames(x$design), collapse = ", "), "\n")
  }
  coefs <- coefficient_names(x)
  if (x$new_design_required && !is.null(coefs)) {
    cat("  new_design columns required:", paste(coefs, collapse = ", "), "\n")
  }
  invisible(x)
}

# Internals --------------------------------------------------------------------

check_marginals <- function(marginals) {
  if (!inherits(marginals, "dependentsimr_marginals")) {
    stop("'marginals' must be the result of fit_marginals()")
  }
}

coefficient_names <- function(marginals) {
  for (ds in marginals$datasets) {
    if (!is.null(ds$coefficients)) return(colnames(ds$coefficients))
  }
  NULL
}

as_dataset_list <- function(data) {
  if (is.matrix(data) || is.data.frame(data)) {
    return(list(datasets = list(data = as_numeric_matrix(data, "data")), single = TRUE))
  }
  if (!is.list(data) || length(data) == 0) {
    stop("'data' must be a matrix, a data frame, or a named list of them")
  }
  if (is.null(names(data)) || any(names(data) == "") || anyDuplicated(names(data))) {
    stop("When 'data' is a list, it must have unique names for each dataset")
  }
  datasets <- lapply(names(data), function(nm) as_numeric_matrix(data[[nm]], nm))
  names(datasets) <- names(data)
  list(datasets = datasets, single = FALSE)
}

as_numeric_matrix <- function(x, name) {
  if (is.data.frame(x)) {
    is_num <- vapply(x, is.numeric, TRUE)
    if (ncol(x) > 1 && !is_num[1] && all(is_num[-1])) {
      # First column holds the feature names
      ids <- as.character(x[[1]])
      if (anyDuplicated(ids) || anyNA(ids)) {
        stop("The first column of '", name, "' is used as feature names but has duplicate or missing values")
      }
      x <- as.matrix(x[, -1, drop = FALSE])
      rownames(x) <- ids
    } else if (all(is_num)) {
      x <- as.matrix(x)
    } else {
      stop(
        "Data frame '", name, "' must be numeric, except optionally its first column of feature names"
      )
    }
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("Dataset '", name, "' must be a numeric matrix or data frame")
  }
  if (anyNA(x)) {
    stop("Dataset '", name, "' contains missing values")
  }
  x
}

check_design <- function(design, n, sample_names = NULL) {
  if (is.null(design)) {
    return(matrix(1, nrow = n, ncol = 1, dimnames = list(NULL, "(Intercept)")))
  }
  if (is.data.frame(design)) {
    stop("'design' must be a numeric matrix, e.g. from model.matrix(), not a data frame")
  }
  if (is.null(dim(design))) {
    stop("'design' must be a matrix with one row per sample")
  }
  if (!is.numeric(design) || anyNA(design)) {
    stop("'design' must be a numeric matrix without missing values")
  }
  if (nrow(design) != n) {
    stop("'design' must have one row per sample (", n, ") but has ", nrow(design))
  }
  if (is.null(colnames(design))) {
    colnames(design) <- paste0("X", seq_len(ncol(design)))
  }
  if (anyDuplicated(colnames(design))) {
    stop("'design' must have unique column names")
  }
  if (!is.null(rownames(design)) && !is.null(sample_names) &&
      !all(rownames(design) == as.character(seq_len(n))) &&
      !identical(rownames(design), sample_names)) {
    stop("The row names of 'design' do not match the sample (column) names of 'data'")
  }
  if (qr(design)$rank < ncol(design)) {
    stop("'design' is not full rank: some columns are linear combinations of others")
  }
  if (ncol(design) >= n) {
    stop("'design' must have fewer columns than there are samples")
  }
  design
}

fit_normal <- function(Y, design) {
  n <- ncol(Y)
  p <- ncol(design)
  fit <- stats::lm.fit(design, t(Y))
  coefficients <- t(matrix(fit$coefficients, nrow = p))
  resid <- t(matrix(fit$residuals, nrow = n))
  sigma2 <- rowSums(resid^2) / (n - p)
  # Residuals scaled to have unit variance; zero for features perfectly fit by the design
  residuals <- resid / sqrt(sigma2)
  residuals[sigma2 <= 0, ] <- 0
  list(coefficients = coefficients, sigma2 = sigma2, residuals = residuals)
}

fit_poisson <- function(Y, design) {
  if (any(Y < 0)) {
    stop("method = 'poisson' requires non-negative data")
  }
  m <- nrow(Y)
  n <- ncol(Y)
  p <- ncol(design)
  coefficients <- matrix(NA_real_, nrow = m, ncol = p)
  mu <- matrix(NA_real_, nrow = m, ncol = n)
  all_zero <- rowSums(Y) == 0

  intercept_only <- p == 1 && all(design[, 1] == design[1, 1])
  if (intercept_only) {
    # Closed form MLE
    lambda <- rowMeans(Y)
    coefficients[!all_zero, 1] <- log(lambda[!all_zero]) / design[1, 1]
    mu[!all_zero, ] <- lambda[!all_zero]
  } else {
    family <- stats::poisson()
    for (i in which(!all_zero)) {
      fit <- suppressWarnings(stats::glm.fit(design, Y[i, ], family = family))
      coefficients[i, ] <- fit$coefficients
      mu[i, ] <- fit$fitted.values
    }
  }
  # Features that are all zero have no fit (NA coefficients) and are simulated as all zero
  residuals <- poisson_rqr(Y, mu)
  residuals[all_zero, ] <- 0
  list(coefficients = coefficients, residuals = residuals)
}

fit_deseq <- function(Y, design) {
  rlang::check_installed("DESeq2", reason = "Using method='DESeq2' requires DESeq2")
  rlang::check_installed("S4Vectors", reason = "Using method='DESeq2' requires S4Vectors")
  rlang::check_installed("SummarizedExperiment", reason = "Using method='DESeq2' requires SummarizedExperiment")

  counts <- as.matrix(Y)
  sample_names <- if (is.null(colnames(counts))) paste0("sample", seq_len(ncol(counts))) else colnames(counts)
  colnames(counts) <- sample_names
  dummy <- data.frame(dummy = rep(1, ncol(counts)), row.names = sample_names)

  # The design is passed as a model matrix ('full') so arbitrary designs work.
  # The formula design is unused and just needs to be valid.
  full <- matrix(design, nrow = nrow(design), dimnames = list(NULL, make.names(colnames(design), unique = TRUE)))
  suppressMessages({
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, dummy, ~1)
    dds <- DESeq2::DESeq(dds, full = full, betaPrior = FALSE, quiet = TRUE)
  })

  # DESeq2 model:
  #   Y_ij ~ NB(mu_ij, alpha_i)
  #   mu_ij = s_j q_ij
  #   log2(q_ij) = x_j beta_i
  # with s_j the size factor of sample j and alpha_i the dispersion of gene i.
  # DESeq2 renames the coefficients, but keeps the order of the design columns
  coefficients <- as.matrix(stats::coef(dds))
  stopifnot(ncol(coefficients) == ncol(design))
  size_factors <- DESeq2::sizeFactors(dds)
  mu <- 2^(coefficients %*% t(design)) * matrix(size_factors, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)

  if (!is.null(SummarizedExperiment::assays(dds)$replaceCounts)) {
    # DESeq2 has replaced outliers, so we use that data
    counts <- SummarizedExperiment::assays(dds)$replaceCounts
  }

  # We use the basic MLE estimates of dispersion from DESeq2 since the moderated
  # dispersion estimates bias the data such that the simulation is not realistic
  # even though over-estimating dispersion (for low read counts) is fine for inference.
  dispersion <- S4Vectors::mcols(dds)$dispGeneEst

  # Discrete distribution, so use randomized quantile residuals
  residuals <- nb_rqr(counts, mu, dispersion)
  # DESeq2 gives NAs where it can't fit (e.g., all zero genes). We replace those with zeros
  residuals[is.na(residuals)] <- 0
  dimnames(residuals) <- NULL

  list(
    coefficients = unname(coefficients),
    dispersion = dispersion,
    size_factors = unname(size_factors),
    residuals = residuals
  )
}

fit_empirical <- function(Y) {
  # Transform by the empirical CDF, with mass 'smeared' out so that the full range of 0-1 is possible
  lower <- t(apply(Y, 1, function(x) ecdf(x)(x - 1e-8)))
  upper <- t(apply(Y, 1, function(x) ecdf(x)(x)))
  p <- matrix(runif(nrow(Y) * ncol(Y), lower, upper), nrow = nrow(Y), ncol = ncol(Y))
  residuals <- qnorm(p)
  dimnames(residuals) <- NULL
  list(coefficients = NULL, values = unname(Y), residuals = residuals)
}
