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

  middle_diag_offsets = seq(from=1, length.out=k-1)
  middle_diagonals = lapply(middle_diag_offsets, function(i) {
    diagonal_end = min(last_col+i, last_row)
    rnorm(n = (diagonal_end - i), sd = sd[i+1:diagonal_end])
  })
  bottom_diag_length = min(num_observations, num_variables-k)
  bottom_diagonal = sqrt(rchisq(n = bottom_diag_length, df = seq(num_variables-k, num_variables-k-bottom_diag_length+1, -1))) * population_sd

  if (num_eigs == 0) {
    # Compute all eigenvalues using dense matrices
    # Build the matrix from its diagonals
    mat <- diag(main_diagonal, nrow=last_row, ncol=last_col)
    for (i in middle_diag_offsets) {
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
      do.call(rbind, lapply(middle_diag_offsets, function(i) offset_diag(i, length(middle_diagonals[[i]])) )),
      offset_diag(k, length(bottom_diagonal))
    )
    vals <- c(main_diagonal, middle_diagonals |> unlist(), bottom_diagonal)
    mat <- Matrix::sparseMatrix(locs[,1], locs[,2], x=vals)
    singular_vals <- sparsesvd::sparsesvd(mat, rank=num_eigs)$d
  }
  return(singular_vals^2)
}


#' Efficiently sample the eigenvalues of a random Wishart matrix with spiked eigenvalues and compute the derivative
#' of the singular values with respect to spiked eigenvalues
#'
#' @param spiked_sd The spiked standard deviations
#' @param num_observations The number of observations (aka samples or columns)
#' @param num_variables The number of variables (aka features or rows)
#' @param population_sd the standard deviation of all non-spiked components (num_variables - length(spiked_sd) of them)
#' @param num_eigs The number of eigenvalues to compute. If 0 compute all of them using dense matrix routines. If greater than zero, use sparse matrices and compute that many top eigenvalues.
#'
#' @returns List with both a vector of random singular values of G where G is a random num_variables x num_observations matrix with iid columns from N(0, Sigma) where Sigma is diagonal with entries spiked_sd and all the remaining are population_sd.
#'  and also the Jacobian, where [[i,j]] is the derivative of the ith singular value with respect to the jth spiked SD
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
#' res = sample_spiked_wishart_and_deriv(
#'   spiked_sd = c(500, 100),
#'   num_variables = 1000,
#'   num_observations = 10-1,
#'   num_eigs = 3
#' )
#' res$singular_vals # singular values (of G, i.e., square roots of the eigenvalues of W = G G^T)
#' res$jacobian # jacobian of the singular values (sqrt of the first eigenvalue) with respect to each of the spiked_sd's
sample_spiked_wishart_and_jac <- function(
    spiked_sd,
    num_observations,
    num_variables,
    population_sd = 1.,
    num_eigs = 0) {

  k = length(spiked_sd)
  if (num_variables<= k) { stop("Must have num_variables larger than k, the number of spiked values") }
  if (num_observations < 1) { stop("num_observations must be at least 1")}
  if (num_eigs == 0) { num_eigs = min(num_observations, num_variables) }

  # The matrix is sparse, (k+1)-diagonal only but also is just completely zero beyond this row/col
  last_row = min(num_variables, num_observations+k)
  last_col = min(num_variables, num_observations)

  # Standard deviations of the the variables - but we only care about the first last_row of them
  # all the others are subsumed into the bottom diagonal's chisquare (and are population_sd)
  sd = c(
    spiked_sd,
    rep(population_sd, last_row - k)
  )

  main_diagonal = sqrt(rchisq(n=last_col, df = seq(last_col, 1, -1)))

  middle_diag_offsets = seq(from=1, length.out=k-1)
  middle_diagonals = lapply(middle_diag_offsets, function(i) {
    diagonal_end = min(last_col+i, last_row)
    rnorm(n = (diagonal_end - i), sd = 1)
  })
  bottom_diag_length = min(num_observations, num_variables-k)
  bottom_diagonal = sqrt(rchisq(n = bottom_diag_length, df = seq(num_variables-k, num_variables-k-bottom_diag_length+1, -1)))

  # Compute the top num_eigs eigenvalues using sparse matrices
  locs <- rbind(
    offset_diag(0, length(main_diagonal)),
    do.call(rbind, lapply(middle_diag_offsets, function(i) offset_diag(i, length(middle_diagonals[[i]])) )),
    offset_diag(k, length(bottom_diagonal))
  )
  vals <- c(main_diagonal, middle_diagonals |> unlist(), bottom_diagonal)
  base_mat <- Matrix::sparseMatrix(locs[,1], locs[,2], x=vals)
  mat <- sd * base_mat
  udv <- sparsesvd::sparsesvd(mat, rank=num_eigs)
  jacobian <- lapply(1:num_eigs, function(i) {
    (udv$u[1:k,i] %*% t(udv$v[1:k,i] ) * base_mat[1:k,1:k]) |> as.matrix() |> rowSums()
  })
  jacobian <- do.call(rbind, jacobian)
  return(list(singular_vals = udv$d, jacobian = jacobian))
}

multi_sample_spiked_wishart <- function(
    count,
    spiked_sd,
    num_observations,
    num_variables,
    population_sd = 1.,
    num_eigs = 0) {
  if (num_eigs == 0) { num_eigs = min(num_observations, num_variables) }
  rank = length(spiked_sd)
  samples <- lapply(
    1:count,
    function(i) { sample_spiked_wishart_and_jac(spiked_sd, num_observations, num_variables, population_sd, num_eigs) }
  )
  eig_means = colMeans(do.call(rbind, lapply(samples, function(x) x$singular_vals)))
  # Extract derivs into a 3D matrix, last index being the iteration
  jac = array(lapply(samples, function(x) x$jacobian) |> unlist(),
              dim = c(num_eigs, rank, count))
  jac_mean = apply(jac, c(1,2), mean) # mean over the iterations
  list(
    singular_vals = eig_means,
    jacobian = jac_mean
  )
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
match_with_spiked_wishart <- function(
    desired_eigenvalues,
    rank,
    num_observations,
    num_variables,
    population_sd=1,
    num_iterations=20,
    num_samples_per_iter=300
    ) {
  num_eigs = length(desired_eigenvalues)
  spiked_sd <- desired_eigenvalues[1:min(rank, num_eigs)] / sqrt(num_observations) # Starting guess, not very accurate
  if (rank > num_eigs) { spiked_sd <- c(spiked_sd, rep(min(spiked_sd), rank - num_eigs))}
  beta <- 0.9
  jac_est <- NULL
  resid_est <- NULL
  for (i in 1:num_iterations) {
    # Iterate towards convergence
    # See what eigenvalues the current values produce
    mean_vals <- multi_sample_spiked_wishart(num_samples_per_iter, spiked_sd, num_observations, num_variables, population_sd, num_eigs = num_eigs)
    eig_means = mean_vals$singular_vals
    jac_mean = mean_vals$jacobian
    if (is.null(jac_est)) {
      jac_est <- jac_mean
    } else {
      jac_est <- beta*jac_mean + (1- beta)*jac_est
    }
    resid <- desired_eigenvalues - eig_means
    if (is.null(resid_est)) {
      resid_est <- resid
    } else {
      resid_est <- beta*resid + (1-beta)*resid_est
    }

    learning_rate = 0.1
    # Levenberg–Marquardt
    lambda = 1
    delta <- solve(t(jac_est) %*% jac_est + lambda*diag(rep(1, rank)), t(jac_est) %*% resid_est) |> t()
    spiked_sd <- spiked_sd + delta
  }
  spiked_sd[1,]
}

