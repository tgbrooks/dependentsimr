set.seed(1)

# Shared example data --------------------------------------------------------
true_Sigma <- matrix(c(
  1,   0.8, 0, 0, 0,   0,
  0.8, 1,   0, 0, 0,   0,
  0,   0,   1, 0, 0,   0,
  0,   0,   0, 1, 0,   0,
  0,   0,   0, 0, 1,   0.3,
  0,   0,   0, 0, 0.3, 1
), nrow = 6, ncol = 6)
norm_data <- t(MASS::mvrnorm(n = 20, mu = rep(0, 6), Sigma = true_Sigma))
pois_data <- qpois(pnorm(norm_data), lambda = c(50, 1000, 2, 10, 0.2, 0.5))

counts <- as.matrix(read_counts[, -1])
rownames(counts) <- read_counts$gene_id
counts <- counts[1:200, ]

group <- factor(rep(c("A", "B"), each = 6))
design <- model.matrix(~ group)

# Basic pipeline ------------------------------------------------------------
test_that("all marginal and dependence methods combine", {
  for (marg in c("normal", "poisson", "empirical")) {
    data <- if (marg == "normal") norm_data else pois_data
    marginals <- fit_marginals(data, method = marg)
    for (dep in c("pca", "spiked Wishart", "corpcor", "independent")) {
      dependence <- fit_dependence(marginals, method = dep, rank = 2)
      sim <- simulate_data(marginals, dependence, n_samples = 15)
      expect_true(is.matrix(sim))
      expect_equal(dim(sim), c(6, 15))
      expect_true(all(is.finite(sim)))
    }
  }
})

test_that("DESeq2 marginals work and keep gene names", {
  marginals <- fit_marginals(counts, method = "DESeq2")
  dependence <- fit_dependence(marginals, method = "pca", rank = 2)
  sim <- simulate_data(marginals, dependence, n_samples = 8, size_factors = rep(c(0.5, 2), 4))
  expect_equal(dim(sim), c(200, 8))
  expect_equal(rownames(sim), rownames(counts))
  expect_true(all(sim >= 0))
  # All-zero genes stay all zero
  zero_genes <- rowSums(counts) == 0
  expect_true(all(sim[zero_genes, ] == 0))
})

test_that("data frames with a feature-name column are accepted", {
  marginals <- fit_marginals(read_counts[1:50, ], method = "poisson")
  expect_equal(marginals$datasets$data$rownames, read_counts$gene_id[1:50])
  expect_error(fit_marginals(data.frame(a = "x", b = "y", c = 1), method = "poisson"), "numeric")
})

test_that("multiple datasets give a list", {
  marginals <- fit_marginals(list(norm = norm_data, pois = pois_data), method = list(norm = "normal", pois = "poisson"))
  sim <- simulate_data(marginals, fit_dependence(marginals, "pca", rank = 2), n_samples = 10)
  expect_named(sim, c("norm", "pois"))
  expect_equal(dim(sim$norm), c(6, 10))
  expect_equal(dim(sim$pois), c(6, 10))
})

test_that("normal residuals with intercept-only design match standardization", {
  marginals <- fit_marginals(norm_data, method = "normal")
  ds <- marginals$datasets$data
  expect_equal(unname(ds$coefficients[, 1]), rowMeans(norm_data))
  expect_equal(ds$sigma2, apply(norm_data, 1, var))
  expect_equal(unname(ds$residuals), unname((norm_data - rowMeans(norm_data)) / apply(norm_data, 1, sd)))
})

test_that("pca factor sizes reproduce the residual principal components", {
  marginals <- fit_marginals(norm_data, method = "normal")
  dep <- fit_dependence(marginals, method = "pca", rank = 2)
  pc_var <- rowSums((dep$cov$u %*% diag(dep$pc_factor_sizes))^2)
  Sigma <- dep$cov$u %*% diag(dep$pc_factor_sizes)^2 %*% t(dep$cov$u) + diag(pmax(dep$var - pc_var, 0))
  expect_equal(diag(Sigma), dep$var)
  expect_equal(drop(t(dep$cov$u[, 1]) %*% Sigma %*% dep$cov$u[, 1]), dep$cov$d[1]^2)
})

test_that("independent dependence and remove_dependence agree", {
  marginals <- fit_marginals(norm_data, method = "normal")
  dep <- remove_dependence(fit_dependence(marginals, "corpcor"))
  expect_equal(dep$method, "independent")
  expect_equal(dim(simulate_data(marginals, dep, n_samples = 4)), c(6, 4))
})

test_that("extreme Poisson values do not produce infinite residuals", {
  # Previously the naive transform gave Inf here
  data <- matrix(c(rep(0, 19), 1000, rep(5, 20)), nrow = 2, byrow = TRUE)
  marginals <- fit_marginals(data, method = "poisson")
  expect_true(all(is.finite(marginals$datasets$data$residuals)))
})

test_that("tail-safe quantiles avoid Inf for large z", {
  z <- matrix(c(-40, -10, 0, 10, 40), 1)
  q <- dependentsimr:::quantile_from_normal(z, qnbinom, list(mu = 100, size = 10))
  expect_true(all(is.finite(q)))
  expect_equal(q[3], qnbinom(0.5, mu = 100, size = 10))
})

# Designs ------------------------------------------------------------------
test_that("design and new_design must be given together", {
  marginals <- fit_marginals(counts, method = "DESeq2", design = design)
  dependence <- fit_dependence(marginals, "pca", rank = 2)
  expect_error(simulate_data(marginals, dependence, n_samples = 4), "new_design")

  plain <- fit_marginals(norm_data, method = "normal")
  expect_error(
    simulate_data(plain, fit_dependence(plain, "pca"), new_design = cbind("(Intercept)" = rep(1, 4))),
    "not given a 'design'"
  )
})

test_that("DESeq2 with a design simulates group differences", {
  # Make the second group have 4x expression in some genes
  shifted <- counts
  de <- 1:20
  shifted[de, 7:12] <- 4 * shifted[de, 7:12]
  marginals <- fit_marginals(shifted, method = "DESeq2", design = design)
  coefs <- marginals$datasets$data$coefficients
  expect_equal(colnames(coefs), c("(Intercept)", "groupB"))
  expressed <- de[rowSums(counts[de, ]) > 100]
  expect_true(all(abs(coefs[expressed, "groupB"] - 2) < 0.5))

  dependence <- fit_dependence(marginals, "pca", rank = 2)
  new_design <- model.matrix(~ group, data.frame(group = factor(rep(c("A", "B"), each = 3))))
  sim <- simulate_data(marginals, dependence, new_design = new_design)
  expect_equal(dim(sim), c(200, 6))
  expect_null(colnames(sim))

  # A design built from differently-named variables is caught by model.matrix's metadata
  other <- model.matrix(~ g, data.frame(g = factor(rep(c("A", "B"), each = 3))))
  colnames(other) <- colnames(design)
  expect_error(simulate_data(marginals, dependence, new_design = other), "constructed differently")
})

test_that("new_design columns are checked", {
  marginals <- fit_marginals(norm_data, method = "normal", design = model.matrix(~ factor(rep(1:2, 10))))
  dependence <- fit_dependence(marginals, "pca", rank = 2)
  expect_error(simulate_data(marginals, dependence, new_design = cbind(a = rep(1, 3), b = 0)), "must match")
  expect_error(simulate_data(marginals, dependence, new_design = cbind(rep(1, 3))), "columns")
  expect_error(
    simulate_data(marginals, dependence, n_samples = 4, new_design = marginals$design),
    "does not match"
  )
})

test_that("design rows must be samples and full rank", {
  expect_error(fit_marginals(norm_data, "normal", design = cbind(1, rep(1, 20))), "full rank")
  expect_error(fit_marginals(norm_data, "normal", design = cbind(1, 1:5)), "one row per sample")
  expect_error(fit_marginals(pois_data, "empirical", design = cbind(rep(1, 20))), "empirical")
})

test_that("a design removes known group structure from the dependence", {
  # Independent features with a strong group effect: without the design, the group
  # effect looks like dependence; with the design, it should be removed.
  set.seed(2)
  n <- 40
  grp <- rep(c(0, 1), each = n / 2)
  m <- 300
  effects <- rnorm(m, sd = 3)
  data <- matrix(rnorm(m * n), m, n) + outer(effects, grp)
  X <- cbind("(Intercept)" = 1, grp = grp)

  dep_plain <- fit_dependence(fit_marginals(data, "normal"), "pca", rank = 1)
  dep_design <- fit_dependence(fit_marginals(data, "normal", design = X), "pca", rank = 1)
  expect_gt(dep_plain$pc_factor_sizes[1], 5)
  # With n - p = 38 iid samples and 300 features, the top sample PC of pure noise is about (1 + sqrt(300/38))
  expect_lt(dep_design$pc_factor_sizes[1], 1 + sqrt(m / (n - 2)) + 1)

  # Residual dependence is still captured: add a shared factor on top of the group effect
  factor_loadings <- c(rep(1, 50), rep(0, m - 50))
  data2 <- data + outer(factor_loadings, rnorm(n, sd = 2))
  dep2 <- fit_dependence(fit_marginals(data2, "normal", design = X), "pca", rank = 1)
  top <- abs(dep2$cov$u[, 1])
  expect_gt(sum(top[1:50]^2), 0.75)
})

test_that("poisson with a design recovers effects", {
  set.seed(3)
  grp <- rep(c(0, 1), each = 10)
  lambda <- outer(rep(c(5, 50), 10), exp(log(3) * grp))
  data <- matrix(rpois(length(lambda), lambda), nrow(lambda))
  data[1, ] <- 0 # all zero feature
  X <- cbind("(Intercept)" = 1, grp = grp)
  marginals <- fit_marginals(data, "poisson", design = X)
  coefs <- marginals$datasets$data$coefficients
  expect_true(all(is.na(coefs[1, ])))
  expect_lt(abs(stats::median(coefs[-1, "grp"]) - log(3)), 0.2)
  sim <- simulate_data(marginals, fit_dependence(marginals, "corpcor"), new_design = X)
  expect_true(all(sim[1, ] == 0))
  expect_true(all(is.finite(sim)))
})

test_that("leverage-one samples give a warning", {
  X <- cbind(1, c(1, rep(0, 19)))
  expect_warning(fit_marginals(norm_data, "normal", design = X), "perfectly fit")
})

test_that("spline designs work when built from the original terms", {
  rlang::check_installed("splines")
  samples <- data.frame(age = seq(20, 80, length.out = 20))
  mf <- model.frame(~ splines::ns(age, 3), samples)
  X <- model.matrix(mf, samples)
  marginals <- fit_marginals(norm_data, "normal", design = X)
  dependence <- fit_dependence(marginals, "pca", rank = 2)
  new_samples <- data.frame(age = c(25, 50, 75))
  new_X <- model.matrix(stats::delete.response(stats::terms(mf)), new_samples)
  sim <- simulate_data(marginals, dependence, new_design = new_X)
  expect_equal(dim(sim), c(6, 3))
})

# add_effect ---------------------------------------------------------------
test_that("add_effect adds DE that requires a new_design", {
  marginals <- fit_marginals(counts, method = "DESeq2")
  dependence <- fit_dependence(marginals, "pca", rank = 2)
  expressed <- rownames(counts)[rowSums(counts) > 1000][1:10]
  marginals_de <- add_effect(marginals, "treatment", setNames(rep(3, 10), expressed))

  expect_error(simulate_data(marginals_de, dependence, n_samples = 4), "add_effect")
  expect_error(add_effect(marginals_de, "treatment", rep(1, 200)), "overwrite")
  expect_error(add_effect(marginals, "treatment", c(nope = 1)), "Unknown feature")
  expect_error(add_effect(marginals, "treatment", rep(1, 5)), "one entry per feature")

  new_design <- cbind("(Intercept)" = 1, treatment = rep(c(0, 1), each = 20))
  sim <- simulate_data(marginals_de, dependence, new_design = new_design)
  lfc <- log2(rowMeans(sim[expressed, 21:40]) / rowMeans(sim[expressed, 1:20]))
  expect_true(all(abs(lfc - 3) < 0.5))
  others <- setdiff(rownames(counts)[rowSums(counts) > 1000], expressed)
  expect_true(all(abs(log2(rowMeans(sim[others, 21:40]) / rowMeans(sim[others, 1:20]))) < 0.5))
})

test_that("add_effect can overwrite a fit coefficient and handles multiple datasets", {
  X <- cbind("(Intercept)" = 1, grp = rep(c(0, 1), each = 10))
  marginals <- fit_marginals(
    list(norm = norm_data, emp = pois_data),
    method = list(norm = "normal", emp = "empirical")
  )
  expect_error(add_effect(marginals, "grp", rep(1, 6)), "'dataset' must be given")
  expect_error(add_effect(marginals, "grp", rep(1, 6), dataset = "emp"), "empirical")
  marginals <- add_effect(marginals, "grp", rep(100, 6), dataset = "norm")
  marginals <- add_effect(marginals, "grp", c(rep(100, 3), rep(0, 3)), dataset = "norm", overwrite = TRUE)
  sim <- simulate_data(marginals, fit_dependence(marginals, "independent"), new_design = X)
  expect_true(all(rowMeans(sim$norm[1:3, 11:20]) > 90))
  expect_true(all(abs(rowMeans(sim$norm[4:6, ])) < 5))
})

test_that("dependence must come from matching marginals", {
  m1 <- fit_marginals(norm_data, "normal")
  m2 <- fit_marginals(counts, "DESeq2")
  expect_error(simulate_data(m1, fit_dependence(m2, "pca"), n_samples = 3), "not fit from these")
})
