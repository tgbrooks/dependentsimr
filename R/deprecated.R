#' Compute structure of dependency from a given data (deprecated)
#'
#' @description
#' **Deprecated.** Use [fit_marginals()] followed by
#' [fit_dependence()] instead, which allows reusing the (slow) marginal fit
#' across dependence methods and supports a design matrix of known covariates.
#'
#' @param datasets A list of data matrices that we will mimic. All datasets have samples in the columns and features (e.g., genes, proteins) in the rows. All datasets must have the same samples in corresponding columns.
#' @param method One of 'pca', 'spiked Wishart', or 'corpcor', for the method of determining a dependency structure
#' @param rank Number of PCA components to approximate the dependence structure of. Only used if method = 'pca'.
#' @param types The marginal distribution types ('normal', 'poisson', 'DESeq2', or 'empirical'), as a list with entries corresponding to datasets. If just a single value is provided, then it is used for all datasets.
#'
#' @returns A random structure element suitable for use with draw_from_multivariate_corr().
#' @export
#' @keywords internal
#'
#' @examples
#' # Simple data to mimic
#' true_means <- c(0, 1, 10, 1000)
#' data <- rpois(4*12, true_means) |> matrix(4, 12)
#'
#' # Simulate draws mimicking that data
#' rs <- get_random_structure(list(data=data), method="pca", rank=2, type="poisson")
#' draws <- draw_from_multivariate_corr(rs, n_samples=30)
get_random_structure <- function(datasets, method, rank=2, types="normal") {
  .Deprecated(msg = paste(
    "get_random_structure() is deprecated.",
    "Use fit_marginals() and then fit_dependence() instead."
  ))
  if (!is.list(datasets) || is.data.frame(datasets)) {
    stop("'datasets' must be a named list of data matrices")
  }
  if (!(method %in% c("pca", "spiked Wishart", "corpcor"))) {
    stop("Parameter 'method' must be one of 'pca', 'spiked Wishart', or 'corpcor'")
  }
  if (length(types) == 1 && !is.list(types)) {
    types <- lapply(datasets, function(x) types[1])
  }
  stopifnot(all(names(datasets) %in% names(types)))
  types <- types[names(datasets)]

  marginals <- fit_marginals(datasets, method = types)
  dependence <- fit_dependence(marginals, method = method, rank = rank)

  # Old-style per-dataset marginal parameters
  old_marginals <- lapply(marginals$datasets, function(ds) {
    switch(ds$method,
      normal = list(var = ds$sigma2, mean = ds$coefficients[, 1]),
      poisson = list(lambda = ifelse(is.na(ds$coefficients[, 1]), 0, exp(ds$coefficients[, 1]))),
      DESeq2 = list(q = 2^ds$coefficients[, 1], dispersion = ds$dispersion),
      empirical = ds$values
    )
  })

  structure <- list(
    type = types,
    cov_method = method,
    marginals = old_marginals,
    n_features = dependence$n_features,
    n_features_list = dependence$n_features_list,
    n_samples = marginals$n_samples,
    transformed_data = lapply(marginals$datasets, function(ds) unname(ds$residuals)),
    rownames = dependence$rownames
  )
  fields <- switch(method,
    pca = c("rank", "cov", "pc_factor_sizes", "var"),
    `spiked Wishart` = c("rank", "cov", "spiked_sd", "population_sd", "var"),
    corpcor = c("indep_part", "dep_part", "var")
  )
  c(structure, unclass(dependence)[fields])
}

#' Draw random samples from the given random structure (deprecated)
#'
#' @description
#' **Deprecated.** Use [simulate_data()] instead.
#'
#' @param random_structure The structure of the random data to draw, computed by get_random_structure()
#' @param n_samples The number of samples to generate
#' @param size_factors Vector of length n_samples of size factors for DESeq2-type data. By default, all samples are given 1. Not used if no datasets has type DESeq2.
#'
#' @return A list of generated datasets. List names correspond to those of the datasets used to generate the random structure. Each dataset has n_samples columns.
#' @export
#' @keywords internal
#'
#' @examples
#' data <- rpois(4*12, c(0, 1, 10, 1000)) |> matrix(4, 12)
#' rs <- get_random_structure(list(data=data), method="pca", rank=2, type="poisson")
#' draws <- draw_from_multivariate_corr(rs, n_samples=30)
draw_from_multivariate_corr <- function(random_structure, n_samples, size_factors=NULL) {
  .Deprecated(msg = paste(
    "draw_from_multivariate_corr() is deprecated.",
    "Use simulate_data() instead."
  ))
  rs <- random_structure
  marginals <- structure(
    list(
      datasets = lapply(names(rs$n_features_list), function(dataname) {
        mg <- rs$marginals[[dataname]]
        ds <- switch(rs$type[[dataname]],
          normal = list(coefficients = cbind(mg$mean), sigma2 = mg$var),
          poisson = list(coefficients = cbind(log(mg$lambda))),
          DESeq2 = list(coefficients = cbind(log2(mg$q)), dispersion = mg$dispersion),
          empirical = list(coefficients = NULL, values = mg),
          stop("Unknown type '", rs$type[[dataname]], "'")
        )
        if (!is.null(ds$coefficients)) colnames(ds$coefficients) <- "(Intercept)"
        ds$method <- rs$type[[dataname]]
        ds$n_features <- rs$n_features_list[[dataname]]
        ds$rownames <- rs$rownames[[dataname]]
        ds
      }) |> stats::setNames(names(rs$n_features_list)),
      design = matrix(1, nrow = rs$n_samples, ncol = 1, dimnames = list(NULL, "(Intercept)")),
      design_supplied = FALSE,
      new_design_required = FALSE,
      single = FALSE,
      n_samples = rs$n_samples
    ),
    class = "dependentsimr_marginals"
  )
  fields <- switch(rs$cov_method,
    pca = c("rank", "cov", "pc_factor_sizes", "var"),
    `spiked Wishart` = c("rank", "cov", "spiked_sd", "population_sd", "var"),
    corpcor = c("indep_part", "dep_part", "var"),
    stop("Unknown method in random_structure")
  )
  dependence <- structure(
    c(
      list(
        method = rs$cov_method,
        n_features = rs$n_features,
        n_features_list = rs$n_features_list,
        rownames = rs$rownames
      ),
      rs[fields]
    ),
    class = "dependentsimr_dependence"
  )
  simulate_data(marginals, dependence, n_samples = n_samples, size_factors = size_factors)
}

#' Remove all dependence
#'
#' Gives a version of the dependence with no dependence between features, so that
#' all simulated features are independent. Equivalent to
#' `fit_dependence(marginals, method = "independent")`.
#'
#' @param random_structure Fit dependence from [fit_dependence()], or (deprecated)
#'   a random structure from [get_random_structure()]
#'
#' @return The dependence (or random structure) with dependency removed, so all data generated from it will be independent.
#' @export
#'
#' @examples
#' data <- rpois(4*12, c(0, 1, 10, 1000)) |> matrix(4, 12)
#' marginals <- fit_marginals(data, method = "poisson")
#' dependence <- fit_dependence(marginals, method = "pca", rank = 2)
#' sim_indep <- simulate_data(marginals, remove_dependence(dependence), n_samples = 30)
remove_dependence <- function(random_structure) {
  if (inherits(random_structure, "dependentsimr_dependence")) {
    keep <- c("n_features", "n_features_list", "rownames", "n_eff")
    return(structure(
      c(list(method = "independent"), unclass(random_structure)[intersect(keep, names(random_structure))]),
      class = "dependentsimr_dependence"
    ))
  }

  new_structure <- random_structure
  if (new_structure$cov_method == "pca") {
    new_structure$pc_factor_sizes <- rep(0, random_structure$rank)
  } else if (new_structure$cov_method == "spiked Wishart") {
    new_structure$spiked_sd <- rep(0, length(new_structure$spiked_sd))
  } else if (new_structure$cov_method == "corpcor") {
    new_structure$dep_part <- 0*new_structure$dep_part
    new_structure$indep_part <- rep(1, length(new_structure$indep_part))
  } else {
    stop("Unknown method in random structure")
  }
  return(new_structure)
}
