# Tests of the deprecated v1 interface
set.seed(0)
## MULTIVARIATE NORMAL DATA ------------------------------------
# Generate the example MV normal data
library(MASS)
true_Sigma = matrix(c(
  1,      0.8,    0,  0,  0,  0,
  0.8,    1,      0,  0,  0,  0,
  0,      0,      1,  0,  0,  0,
  0,      0,      0,  1,  0,  0,
  0,      0,      0,  0,  1,  0.3,
  0,      0,      0,  0,  0.3, 1
), nrow=6, ncol=6)
norm_data <- t(mvrnorm(n=20, mu=c(0,0,0,0,0,0), Sigma=true_Sigma))

# Simulate draws mimicking that data
rs_normal <- suppressWarnings(get_random_structure(list(data=norm_data), method="pca", rank=2, type="normal"))
draws_normal <- suppressWarnings(draw_from_multivariate_corr(rs_normal, n_samples=30))

test_that("normal draws work", {
  expect_equal(dim(draws_normal$data), c(6,30))
})

# Compare the covariance matrix used to simulate with the expected properties
# that it matches the desired diagonals and principal components
pc_var <- apply((rs_normal$cov$u %*% diag(rs_normal$pc_factor_sizes))^2, 1, sum) # variance we already accounted for
missing_var <- pmax(rs_normal$var - pc_var, 0) # variance left to include as purely independent
Sigma <- rs_normal$cov$u %*% diag(rs_normal$pc_factor_sizes)^2 %*% t(rs_normal$cov$u) + diag(missing_var)
Var <- rs_normal$transformed_data$data %*% t(rs_normal$transformed_data$data) / (ncol(rs_normal$transformed_data$data) - 1)
test_that("pc size factors are accurate", {
  expect_equal(diag(Sigma), rs_normal$var)
  expect_equal((t(rs_normal$cov$u[,1]) %*% Sigma %*% rs_normal$cov$u[,1])[1,],
               rs_normal$cov$d[1]^2)
  expect_equal((t(rs_normal$cov$u[,2]) %*% Sigma %*% rs_normal$cov$u[,2])[1,],
               rs_normal$cov$d[2]^2)
})

## POISSON DATA ------------------------------------------------
# Generate dependent Poisson data from the MV normal data
pois_data <- qpois(pnorm(norm_data), lambda = c(50, 1000, 2, 10, 0.2, 0.5))

# Simulate draws mimicking that data
rs_poisson <- suppressWarnings(get_random_structure(list(data=pois_data), method="pca", rank=1, type="poisson"))
draws_poisson <- suppressWarnings(draw_from_multivariate_corr(rs_poisson, n_samples=30))

test_that("poisson draws work", {
  expect_equal(dim(draws_poisson$data), c(6,30))
  expect(all(draws_poisson$data >= 0), "draws_poisson$data must be all non-negative")
})

## DESEQ2 DATA -------------------------------------------------
# A small amount of RNA-seq read count data
d <- c(
  3563,      3839,      3407,      3745,      4063,      3089,      3219,      4577,      3769,
  0,         0,         0,         0,         0,         0,         0,         0,         0,
  0,         0,         0,         0,         0,         0,         0,         0,         0,
  3626,      2630,      5470,      4230,      4739,      4499,      5862,      6524,      3621,
  317,       310,       254,       449,       431,       463,       292,       325,       299,
  0,         0,         0,         0,         0,         0,         0,         0,         0,
  0,         0,         0,         0,         0,         0,         0,         0,         0,
  0,         0,         0,         0,         0,         0,         0,         0,         0,
  0,         0,         0,         0,         0,         0,         0,         0,         0,
  3,         2,         0,         1,         5,         3,         9,         1,         2,
  2,         1,         0,         2,         0,         0,         0,         0,         0,
  0,         0,         0,         0,         0,         0,         0,         0,         0,
  1,         0,         0,         0,         1,         0,         0,         0,         0,
  1,         1,         2,         0,         0,         1,         0,         0,         1,
  0,         0,         0,         0,         0,         0,         0,         0,         0,
  984,      1014,       777,      1211,      1174,      1096,       973,      1193,      1083,
  0,         0,         0,         0,         0,         0,         0,         0,         0,
  16,        15,        16,         9,        23,        10,        25,        10,         8,
  4473,      3954,      4476,      3965,      4606,      3788,      4084,      4316,      3658,
  0,         0,         0,         0,         0,         0,         0,         1,         0
)
read_counts <- matrix(d, nrow=20, ncol=9, byrow=TRUE)

# Simulate draws mimicking that data
rs_deseq <- suppressWarnings(get_random_structure(list(data=read_counts), method="pca", rank=2, type="DESeq2"))
draws_deseq <- suppressWarnings(draw_from_multivariate_corr(rs_deseq, n_samples=30))

test_that("DESeq2 draws work", {
  expect_equal(dim(draws_deseq$data), c(20,30))
  expect(all(draws_deseq$data >= 0), "draws_deseq$data must be all non-negative")
})


## Empirical fit --------------------------------------------------------
# Using the same as the DESeq2 data
rs_emp <- suppressWarnings(get_random_structure(list(data=read_counts), method="pca", rank=2, type="empirical"))
draws_emp <- suppressWarnings(draw_from_multivariate_corr(rs_emp, n_samples=30))

test_that("empirical draws work", {
  expect_equal(dim(draws_emp$data), c(20,30))
  expect_contains(read_counts, draws_emp$data)
})

## Multi-omics fit ------------------------------------------------------
# Example using more than one marginal data type, such as might occur when
# simulating multi-omics where each omics mode has a different type
# of marginal distribution.
datasets <- list(
  norm = norm_data,
  pois = pois_data
)
types <- list(
  norm = "normal",
  pois = "poisson"
)

rs_multi <- suppressWarnings(get_random_structure(datasets, method="pca", rank=2, types=types))
draws_multi <- suppressWarnings(draw_from_multivariate_corr(rs_multi, n_samples=10))

test_that("multiple datasets works", {
  expect_equal(dim(draws_multi$norm), c(6,10))
  expect_equal(dim(draws_multi$pois), c(6,10))
})

test_that("old interface gives deprecation warnings", {
  expect_warning(get_random_structure(list(data=norm_data), method="pca", rank=2, type="normal"), "deprecated")
  expect_warning(draw_from_multivariate_corr(rs_normal, n_samples=3), "deprecated")
})

test_that("modifying old-style marginals still affects draws", {
  rs_mod <- rs_deseq
  rs_mod$marginals$data$q[1] <- 1e6 * rs_mod$marginals$data$q[1]
  draws <- suppressWarnings(draw_from_multivariate_corr(rs_mod, n_samples=10))$data
  expect_true(all(draws[1,] > 1e8))
})

test_that("remove_dependence works on old structures", {
  rs_indep <- remove_dependence(rs_normal)
  expect_equal(rs_indep$pc_factor_sizes, c(0, 0))
  draws <- suppressWarnings(draw_from_multivariate_corr(rs_indep, n_samples=5))
  expect_equal(dim(draws$data), c(6, 5))
})

test_that("remove_dependence on an old spiked Wishart structure gives independent draws", {
  rs <- suppressWarnings(get_random_structure(list(data=norm_data), method="spiked Wishart", rank=2, type="normal"))
  rs_indep <- remove_dependence(rs)
  dep <- structure(
    list(method = "spiked Wishart", n_features = 6, rank = 2, cov = rs_indep$cov,
         spiked_sd = rs_indep$spiked_sd, population_sd = rs_indep$population_sd),
    class = "dependentsimr_dependence"
  )
  N <- 20000
  C <- tcrossprod(dependentsimr:::draw_latent(dep, N)) / N
  expect_equal(C, diag(rs$population_sd^2, 6), tolerance = 0.05)
})
