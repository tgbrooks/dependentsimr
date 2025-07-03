offset_diag <- function(offset, length) {
  cbind(
    seq(offset+1, offset+length),
    seq(1, length)
  )
}

#' Efficiently sample the singular values corresponding to a random Wishart matrix with spiked eigenvalues
#' Specifically, if W = G G^T with each column of G drawn iid from N(0, Sigma), then W is a Wishart matrix and
#' this function samples the singular values of G. The eigenvalues of W are just the squares of the singular values.
#' Here, Sigma is diagonal with its leading entries from spiked_sd^2 and all remaining entries are population_sd^2.
#'
#' @param spiked_sd The spiked standard deviations
#' @param num_observations The number of observations (aka samples or columns)
#' @param num_variables The number of variables (aka features or rows)
#' @param population_sd the standard deviation of all non-spiked components (num_variables - length(spiked_sd) of them)
#' @param num_eigs The number of eigenvalues to compute. If 0 compute all of them using dense matrix routines. If greater than zero, use sparse matrices and compute that many top eigenvalues.
#'
#' @returns Vector of random singular values of G where G is a random num_variables x num_observations matrix with iid columns from N(0, Sigma) where Sigma is diagonal with entries spiked_sd^2 and all the remaining are population_sd^2.
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
  return(singular_vals)
}

# Compute svd using sparsesvd but return a guaranteed number of singular values
# since sparsesvd will truncate whenever an svd happens to small even if we set tol=0
sparsesvd <- function(mat, rank) {
  udv <- sparsesvd::sparsesvd(mat, rank=rank)

  if (length(udv$d) < rank) {
    missing <- rank - length(udv$d)
    num_vars <- dim(udv$u)[[1]]
    num_obs <- dim(udv$v)[[1]]
    udv$d <- c(udv$d, rep(0, missing))
    udv$u <- cbind(udv$u, rep(0, missing*num_vars) |> matrix(num_vars, missing))
    udv$v <- cbind(udv$v, rep(0, missing*num_obs) |> matrix(num_obs, missing))
  }
  udv
}



#' Efficiently sample the singular values corresponding to a random Wishart matrix with spiked eigenvalues and the Jacobian
#' I.e., these are the singular values of G if GG^T is Wishart. The square of these give the eigenvalues of the random Wishart matrix.
#'
#' @param spiked_sd The spiked standard deviations
#' @param num_observations The number of observations (aka samples or columns)
#' @param num_variables The number of variables (aka features or rows)
#' @param population_sd the standard deviation of all non-spiked components (num_variables - length(spiked_sd) of them)
#' @param num_eigs The number of eigenvalues to compute. If 0 compute all of them using dense matrix routines. If greater than zero, use sparse matrices and compute that many top eigenvalues.
#'
#' @returns List with a vector of random singular values of G where G is a random num_variables x num_observations matrix with iid columns from N(0, Sigma) where Sigma is diagonal with entries spiked_sd^2 and all the remaining are population_sd^2.
#'  and also the Jacobian, where `\[i,j\]` is the derivative of the ith singular value with respect to the jth spiked SD and the
#'  gradient of the population_sd variable
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
#' res = sample_spiked_wishart_and_jac(
#'   spiked_sd = c(500, 100),
#'   num_variables = 1000,
#'   num_observations = 10-1,
#'   num_eigs = 3
#' )
#' res$singular_vals # singular values of G, (i.e., square roots of the eigenvalues of W = G G^T)
#' res$jacobian # jacobian of the singular values with respect to the spiked_sd's
#' res$pop_sd_grad # gradient of population_sd parameter
sample_spiked_wishart_and_jac <- function(
    spiked_sd,
    num_observations,
    num_variables,
    population_sd = 1.,
    num_eigs = 0) {

  k = length(spiked_sd)
  if (num_variables <= k) { stop("Must have num_variables larger than k, the number of spiked values") }
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
  udv <- sparsesvd(mat, rank=num_eigs)
  jacobian <- lapply(1:num_eigs, function(i) {
    udv$u[1:k,i] * Matrix::rowSums(base_mat[1:k,] %*% udv$v[,i, drop=FALSE])
  })
  jacobian <- do.call(rbind, jacobian)

  pop_indexes <- (k+1):nrow(udv$u) # the indexes corresponding to population_sd
  pop_sd_grad <- lapply(1:num_eigs, function(i) {
    (udv$u[pop_indexes,i] * Matrix::rowSums(base_mat[pop_indexes,] %*% udv$v[,i,drop=FALSE])) |> sum()
  }) |> unlist()

  return(list(singular_vals = udv$d, jacobian = jacobian, pop_sd_grad = pop_sd_grad))
}

#' Compute means of each singular value and the mean Jacobian, see sample_spiked_wishart_and_jac
#'
#' @param count The number of samples to compute the mean of
#' @param spiked_sd The spiked standard deviations
#' @param num_observations The number of observations (aka samples or columns)
#' @param num_variables The number of variables (aka features or rows)
#' @param population_sd the standard deviation of all non-spiked components (num_variables - length(spiked_sd) of them)
#' @param num_eigs The number of eigenvalues to compute. If 0 compute all of them using dense matrix routines. If greater than zero, use sparse matrices and compute that many top eigenvalues.
#'
#' @returns List with a vector of mean singular values of G where G is a random num_variables x num_observations matrix with iid columns from N(0, Sigma) where Sigma is diagonal with entries spiked_sd^2 and all the remaining are population_sd^2.
#'  and also the mean Jacobian, where `\[i,j\]` is the derivative of the ith singular value with respect to the jth spiked SD, and the
#'  gradient of the population_sd parameter
#'
#' @importFrom stats rnorm
#' @importFrom stats rchisq
#' @export
#'
#' @examples
#' # Sample 10 times from the spiked Wishart distribution with (500, 100, 1, ..., 1) singular values
#' # and take the means of the singular values as well as derivatives (jacobian and pop_sd_grad)
#' mean_vals <- multi_sample_spiked_wishart(
#'     count = 10,
#'     spiked_sd = c(500, 100),
#'     num_observations = 10-1,
#'     num_variables = 1000,
#'     num_eigs = 3
#' )
#' mean_vals$singular_vals
#' mean_vals$jacobian
#' mean_vals$pop_sd_grad
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

  # Extract derivs into a 3D array, last index being the iteration
  jac = array(lapply(samples, function(x) x$jacobian) |> unlist(),
              dim = c(num_eigs, rank, count))
  jac_mean = apply(jac, c(1,2), mean) # mean over the iterations
  # Similarly for population SD's gradient
  pop_sd_grad = array(lapply(samples, function(x) x$pop_sd_grad) |> unlist(), dim= c(num_eigs, count))
  pop_sd_grad_mean <- apply(pop_sd_grad, 1, mean)

  list(
    singular_vals = eig_means,
    jacobian = jac_mean,
    pop_sd_grad = pop_sd_grad_mean
  )
}

#' Compute what spiked SD values will give you the desired top eigenvalues by iteratively solving
#'
#' @param desired_eigenvalues Vector of top eigenvalues of sample covariance matrix that should be approximated
#' @param rank number of 'spiked' dimensions to fit
#' @param num_observations Number of observations (samples, columns) in the data matrix
#' @param num_variables Number of variables (features, rows) in the data matrix
#' @param population_sd Standard deviation for allnon-spiked dimensions
#' @param num_iterations Number of iterations to perform to to optimize. Increase to get a closer fit
#' @param num_samples_per_iter Number of eigenvalues to samples to estimate mean eigenvalues for each iteration.
#'    Increase to get a closer fit
#'
#' @export
#' @returns Fit values of spiked SDs that approximately match the desired eigenvalues
#'
#' @examples
#' # First, simulate a dataset with known true values to see if we can get close
#' true_spiked_sd <- c(500, 100)
#' num_variables <- 1000
#' num_observations <- 10-1
#' sampled_eigenvalues <- sample_spiked_wishart(true_spiked_sd, num_observations,
#'       num_variables, num_eigs=2)
#' # Fit some spiked SDs that give eigenvalues similar to sampled_eigenvalues
#' fit <- match_with_spiked_wishart(sampled_eigenvalues, rank=2, num_observations, num_variables,
#'           num_iterations=5, num_samples_per_iter=50)
#' # fit$spiked_sd should now be close to giving the specified sampled_eigenvalues (in expectation)
#' # fit$spiked_sd won't match the original true_spiked_sd too closely since it only
#' # fits the single sample # that gave  sampled_eigenvalue
#' # fit$population_sd gives fit value for the population SD
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
    pop_grad = mean_vals$pop_sd_grad
    if (is.null(jac_est)) {
      jac_est <- jac_mean
      pop_grad_est <- pop_grad
    } else {
      jac_est <- beta*jac_mean + (1-beta)*jac_est
      pop_grad_est <- beta*pop_grad + (1-beta)*pop_grad_est
    }
    resid <- desired_eigenvalues - eig_means
    if (is.null(resid_est)) {
      resid_est <- resid
    } else {
      resid_est <- beta*resid + (1-beta)*resid_est
    }

    # combined Jacobian for derivatives WRT the spiked SDs and the population SD
    combined <- cbind(jac_est, pop_grad_est)

    # Levenberg–Marquardt
    lambda = 1
    delta <- solve(t(combined) %*% combined + lambda*diag(rep(1, rank+1)), t(combined) %*% resid_est) |> t()
    spiked_sd <- spiked_sd + delta[1:rank]
    population_sd <- population_sd + delta[rank+1]
  }
  list(spiked_sd=spiked_sd, population_sd=population_sd)
}

