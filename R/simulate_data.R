#' Simulate data from fit marginals and dependence
#'
#' The final step of a simulation: draw new samples whose features have the
#' marginal distributions from [fit_marginals()] and the dependence between
#' features from [fit_dependence()]. Samples are independent of each other.
#'
#' If a `design` was given to [fit_marginals()] (or coefficients were added with
#' [add_effect()]), then `new_design` is required and gives the covariates of each
#' simulated sample. It must have exactly the same columns, in the same order, as
#' the coefficients (i.e., the original design plus any added effects) and must be
#' constructed in the same way as the original design. In particular, for
#' covariates with data-dependent encodings such as splines (`splines::ns()`) or
#' polynomials (`poly()`), build `new_design` from the terms of the original model
#' frame, e.g. `model.matrix(delete.response(terms(model.frame(~ ns(age, 3), samples))), new_samples)`,
#' rather than calling `model.matrix(~ ns(age, 3), new_samples)` directly, which
#' would silently choose different knots.
#'
#' @param marginals Fit marginals from [fit_marginals()], possibly modified by [add_effect()]
#' @param dependence Fit dependence from [fit_dependence()] on the same marginals
#' @param n_samples The number of samples to generate. Defaults to `nrow(new_design)` if that is given.
#' @param new_design Design matrix for the simulated samples, one row per sample.
#'   Required if and only if a `design` was given to [fit_marginals()] or effects
#'   were added with [add_effect()].
#' @param size_factors Vector of length `n_samples` of size factors (relative
#'   library sizes) for DESeq2-method data. Defaults to all 1, which corresponds to
#'   the average read depth of the original data. Ignored for other methods.
#'
#' @return The simulated data, with features in rows and `n_samples` samples in
#'   columns. A matrix if `fit_marginals()` was given a single matrix or data frame,
#'   otherwise a list of matrices with the same names as the original datasets.
#' @export
#'
#' @examples
#' # Generate example data with Sigma as its covariance matrix
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
#' marginals <- fit_marginals(norm_data, method = "normal")
#' dependence <- fit_dependence(marginals, method = "pca", rank = 2)
#' sim <- simulate_data(marginals, dependence, n_samples = 30)
#'
#' # With a design: data from two groups
#' group <- factor(rep(c("A", "B"), each = 10))
#' design <- model.matrix(~ group)
#' marginals <- fit_marginals(norm_data, method = "normal", design = design)
#' dependence <- fit_dependence(marginals, method = "pca", rank = 2)
#' new_design <- model.matrix(~ group, data.frame(group = factor(rep(c("A", "B"), each = 15))))
#' sim <- simulate_data(marginals, dependence, new_design = new_design)
simulate_data <- function(marginals, dependence, n_samples = NULL, new_design = NULL, size_factors = NULL) {
  check_marginals(marginals)
  check_dependence(dependence)
  if (!identical(unname(marginals_fingerprint(marginals)), unname(dependence_fingerprint(dependence)))) {
    stop("'dependence' was not fit from these 'marginals': the datasets or features differ")
  }

  coef_names <- coefficient_names(marginals)
  if (is.null(coef_names)) coef_names <- colnames(marginals$design) # all empirical

  # Determine the design of the new samples
  if (is.null(new_design)) {
    if (marginals$new_design_required) {
      if (marginals$design_supplied) {
        stop("'new_design' is required since a 'design' was given to fit_marginals()")
      } else {
        stop("'new_design' is required since effects were added with add_effect(); it needs columns: ",
             paste(coef_names, collapse = ", "))
      }
    }
    if (is.null(n_samples)) {
      stop("Either 'n_samples' or 'new_design' must be given")
    }
    new_design <- matrix(1, nrow = n_samples, ncol = 1, dimnames = list(NULL, coef_names))
  } else {
    if (!marginals$new_design_required) {
      stop(
        "'new_design' was given but fit_marginals() was not given a 'design'. ",
        "Provide a 'design' to fit_marginals() or use add_effect() to add coefficients"
      )
    }
    new_design <- check_new_design(new_design, marginals, coef_names)
    if (is.null(n_samples)) {
      n_samples <- nrow(new_design)
    } else if (n_samples != nrow(new_design)) {
      stop("'n_samples' (", n_samples, ") does not match the number of rows of 'new_design' (", nrow(new_design), ")")
    }
  }
  if (!(is.numeric(n_samples) && length(n_samples) == 1 && n_samples >= 1 && n_samples == round(n_samples))) {
    stop("'n_samples' must be a positive integer")
  }

  if (is.null(size_factors)) {
    size_factors <- rep(1, n_samples)
  } else if (length(size_factors) != n_samples) {
    stop("'size_factors' must have length n_samples (", n_samples, ")")
  }

  # Name the simulated samples by the rows of new_design, unless those are just model.matrix()'s 1, 2, 3, ...
  sample_names <- rownames(new_design)
  if (identical(sample_names, as.character(seq_len(n_samples)))) sample_names <- NULL

  latent <- draw_latent(dependence, n_samples)

  start_row <- 1
  draws <- list()
  for (dataname in names(marginals$datasets)) {
    ds <- marginals$datasets[[dataname]]
    m <- ds$n_features
    end_row <- start_row + m - 1
    z <- latent[start_row:end_row, , drop = FALSE]
    start_row <- end_row + 1

    if (ds$method != "empirical") {
      linear <- ds$coefficients %*% t(new_design) # features x samples
    }
    # Adjust the normal draws to have the correct marginals
    if (ds$method == "normal") {
      d <- linear + sqrt(ds$sigma2) * z
    } else if (ds$method == "poisson") {
      d <- quantile_from_normal(z, stats::qpois, list(lambda = exp(linear)))
      d[is.na(d)] <- 0 # NAs come from all-zero features in the original data
    } else if (ds$method == "DESeq2") {
      s <- matrix(size_factors, nrow = m, ncol = n_samples, byrow = TRUE)
      mu <- s * 2^linear
      d <- quantile_from_normal(z, stats::qnbinom, list(mu = mu, size = 1 / ds$dispersion))
      d[is.na(d)] <- 0 # NAs come from all-zero features in the original data
    } else if (ds$method == "empirical") {
      p <- pnorm(z)
      d <- matrix(NA_real_, nrow = m, ncol = n_samples)
      for (i in seq_len(m)) {
        d[i, ] <- quantile(ds$values[i, ], p[i, ], type = 1, names = FALSE)
      }
    } else {
      stop("Unknown marginal method '", ds$method, "'")
    }
    d <- matrix(d, nrow = m, ncol = n_samples)
    rownames(d) <- ds$rownames
    colnames(d) <- sample_names
    draws[[dataname]] <- d
  }

  if (marginals$single) draws[[1]] else draws
}

check_new_design <- function(new_design, marginals, coef_names) {
  if (is.data.frame(new_design)) {
    stop("'new_design' must be a numeric matrix, e.g. from model.matrix(), not a data frame")
  }
  if (is.null(dim(new_design)) || !is.numeric(new_design) || anyNA(new_design)) {
    stop("'new_design' must be a numeric matrix without missing values")
  }
  if (is.null(colnames(new_design))) {
    if (ncol(new_design) != length(coef_names)) {
      stop("'new_design' must have ", length(coef_names), " columns: ", paste(coef_names, collapse = ", "))
    }
    colnames(new_design) <- coef_names
  }
  if (!identical(colnames(new_design), coef_names)) {
    stop(
      "The columns of 'new_design' must match the fit coefficients, in order.\n",
      "  Expected: ", paste(coef_names, collapse = ", "), "\n",
      "  Got:      ", paste(colnames(new_design), collapse = ", ")
    )
  }
  # model.matrix() records how columns map to terms and the contrasts used for factors.
  # If both designs came from model.matrix(), these must agree.
  design <- marginals$design
  if (ncol(design) == ncol(new_design)) {
    for (a in c("assign", "contrasts")) {
      old <- attr(design, a)
      new <- attr(new_design, a)
      if (!is.null(old) && !is.null(new) && !identical(old, new)) {
        stop("'new_design' was constructed differently than 'design' (its '", a, "' attribute differs)")
      }
    }
  }
  new_design
}

marginals_fingerprint <- function(marginals) {
  list(
    lapply(marginals$datasets, function(ds) ds$n_features),
    lapply(marginals$datasets, function(ds) ds$rownames)
  )
}

dependence_fingerprint <- function(dependence) {
  list(dependence$n_features_list, dependence$rownames)
}
