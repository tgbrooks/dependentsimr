offset_diag <- function(offset, length) {
  cbind(
    seq(offset+1, offset+length),
    seq(1, length)
  )
}

#' Efficiently sample the eigenvalues of a random Wishart matrix with spiked eigenvalues
#'
#' @param spiked_sd The spiked standard deviations
#' @param num_observations The number of observations (aka samples or columns)
#' @param num_variables The number of variables (aka features or rows)
#' @param population_sd the standard deviation of all non-spiked components (num_variables - length(spiked_sd) of them)
#' @param num_eigs The number of eigenvalues to compute. If 0 compute all of them using dense matrix routines. If greater than zero, use sparse matrices and compute that many top eigenvalues.
#'
#' @returns Vector of random eigenvalues of W = G G^T where G is a random num_variables x num_observations matrix with iid columns from N(0, Sigma) where Sigma is diagonal with entries spiked_sd and all the remaining are population_sd.
#'
#' @importFrom stats rnorm
#' @importFrom stats rchisq
#' @export
#'
#' @examples
#' set.seed(0)
#' # Sample eigenvalues of a covariance matrix of a 10 sample study with 1000 variables such that
#' # the top two (underlying true distribution of the data, not sample) principal components
#' # have SDs of 100 and 10 and the remaining 98 have 1
#' sample_spiked_wishart(
#'   spiked_sd = c(500, 100),
#'   num_variables = 1000,
#'   num_observations = 10-1,
#'   num_eigs = 3
#' )
#'
sample_spiked_wishart <- function(
    spiked_sd,
    num_observations,
    num_variables,
    population_sd = 1.,
    num_eigs = 0) {

  k = length(spiked_sd)
  if (num_variables<= k) { stop("Must have num_variables larger than k, the number of spiked values") }
  if (num_observations < 1) { stop("num_observations must be at least 1")}

  # The matrix is sparse, (k+1)-diagonal only but also is just completely zero beyond this row/col
  last_row = min(num_variables, num_observations+k)
  last_col = min(num_variables, num_observations)

  # Standard deviations of the the variables - but we only care about the first last_row of them
  # all the others are subsumed into the bottom diagonal's chisquare (and are population_sd)
  sd = c(
    spiked_sd,
    rep(population_sd, last_row - k)
  )

  main_diagonal = sd[1:last_col] * sqrt(rchisq(n=last_col, df = seq(last_col, 1, -1)))

  middle_diagonals = lapply(seq(2,k), function(i) {
    diagonal_end = min(last_col+i-1, last_row)
    rnorm(n = (diagonal_end - i + 1), sd = sd[i:diagonal_end])
  })
  bottom_diag_length = min(num_observations, num_variables-k)
  bottom_diagonal = sqrt(rchisq(n = bottom_diag_length, df = seq(num_variables-k, num_variables-k-bottom_diag_length+1, -1))) * population_sd

  if (num_eigs == 0) {
    # Compute all eigenvalues using dense matrices
    # Build the matrix from its diagonals
    mat <- diag(main_diagonal, nrow=last_row, ncol=last_col)
    for (i in 1:(k-1)) {
      locs <- offset_diag(i, length(middle_diagonals[[i]]))
      mat[locs] <- middle_diagonals[[i]]
    }
    locs <- offset_diag(k,length(bottom_diagonal))
    mat[locs] <- bottom_diagonal

    singular_vals <- svd(mat, nu=0, nv=0)$d
  } else {
    # Compute the top num_eigs eigenvalues using sparse matrices
    locs <- rbind(
      offset_diag(0, length(main_diagonal)),
      do.call(rbind, lapply(1:(k-1), function(i) offset_diag(i, length(middle_diagonals[[i]])) )),
      offset_diag(k, length(bottom_diagonal))
    )
    vals <- c(main_diagonal, middle_diagonals |> unlist(), bottom_diagonal)
    mat <- Matrix::sparseMatrix(locs[,1], locs[,2], x=vals)
    singular_vals <- sparsesvd::sparsesvd(mat, rank=num_eigs)$d
  }
  return(singular_vals^2)
}

#' Compute what spiked SD values will give you the desired top eigenvalues by iteratively solving
#'
#' @param desired_eigenvalues Vector of top eigenvalues of sample covariance matrix that should be approximated
#' @param num_observations Number of observations (samples, columns) in the data matrix
#' @param num_variables Number of variables (features, rows) in the data matrix
#' @param num_iterations Number of iterations to perform to to optimize. Increase to get a closer fit
#' @param num_samples_per_iter Number of eigenvalues to samples to estimate mean eigenvalues for each iteration.
#'    Increase to get a closer fit
#'
#' @export
#' @returns Fit values of spiked SDs that approximately match the desired eigenvalues
#'
#' @examples
#' # First, simulate a dataset with known true values to see if we can get close
#' true_spiked_sd = c(500, 100)
#' num_variables = 1000
#' num_observations = 10-1
#' sampled_eigenvalues <- sample_spiked_wishart(true_spiked_sd, num_variables,
#'       num_observations, num_eigs=2)
#' # Fit some spiked SDs that give eigenvalues similar to sampled_eigenvalues
#' fit_spiked_sd <- match_with_spiked_wishart(sampled_eigenvalues, num_observations, num_variables)
#' # fit spiked_sd should now be close to giving the specified sampled_eigenvalues (in expectation)
#' # Fit spiked_sd won't match the original true_spiked_sd too closely since it only
#' # fits the single sample # that gave  sampled_eigenvalues
match_with_spiked_wishart <- function(desired_eigenvalues, num_observations, num_variables, num_iterations=10, num_samples_per_iter=300) {

  spiked_sd <- sqrt(desired_eigenvalues) # Starting guess, not very accurate
  for (i in 1:num_iterations) {
    # Iterate towards convergence
    # See what eigenvalues the current values produce
    samples <- replicate(num_samples_per_iter, sample_spiked_wishart(spiked_sd, num_observations, num_variables, num_eigs = length(desired_eigenvalues)))
    eig_means = rowMeans(samples)

    ratios <- desired_eigenvalues / eig_means
    exponent <- 1/sqrt(i) # Dampens the updates as we go further in iterations to converge better

    # Update towards the target desired eigenvalues
    spiked_sd <- sort(spiked_sd * ratios^exponent, decreasing=TRUE)

  }
  spiked_sd
}

