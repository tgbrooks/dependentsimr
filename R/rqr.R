# Randomized quantile residuals and tail-safe quantile functions

#' @keywords internal
log_add <- function(a, b) {
    # log(exp(a) + exp(b)), vectorized, safe when one side is -Inf
    m <- pmax(a, b)
    d <- -abs(a - b)
    d[is.na(d)] <- 0 # both -Inf
    out <- m + log1p(exp(d))
    out[is.infinite(m) & m < 0] <- -Inf
    out
}

#' @keywords internal
#' @importFrom stats pnbinom
#' @importFrom stats dnbinom
#' @importFrom stats runif
#' @importFrom stats qnorm
nb_rqr <- function(counts, mu, disp) {
    # randomized quantile residuals for negative binomial (DESeq2 model)
    # The naive approach sometimes yields infinites that we can't work with
    # this splits by low and high tails and works in log space
    y <- as.vector(counts)
    m <- as.vector(mu)
    size <- 1 / as.vector(disp)       # recycles if disp is per-row/col vector of right length
    if (length(size) != length(y)) {
        size <- rep_len(size, length(y))
    }
    u <- runif(length(y))

    z <- numeric(length(y))
    bad <- is.na(m) | is.na(size) | is.na(y)
    low <- !bad & (y <= m) # lower-tail side
    upp <- !bad & (y >  m) # upper-tail side

    # lower side: p = F(y-1) + U * P(Y=y), z = qnorm(p)
    if (any(low)) {
        yl <- y[low]; ml <- m[low]; sl <- size[low]
        lF  <- pnbinom(yl - 1, mu = ml, size = sl, lower.tail = TRUE, log.p = TRUE)
        ld  <- dnbinom(yl,     mu = ml, size = sl, log = TRUE)
        lp  <- log_add(lF, log(u[low]) + ld)
        z[low] <- qnorm(lp, log.p = TRUE)
    }

    # upper side: s = S(y) + U * P(Y=y), z = -qnorm(s)
    if (any(upp)) {
        yu <- y[upp]; mu_ <- m[upp]; su <- size[upp]
        lS  <- pnbinom(yu, mu = mu_, size = su, lower.tail = FALSE, log.p = TRUE)
        ld  <- dnbinom(yu, mu = mu_, size = su, log = TRUE)
        ls  <- log_add(lS, log(u[upp]) + ld)
        z[upp] <- qnorm(ls, lower.tail = FALSE, log.p = TRUE)
    }

    z[bad] <- 0
    dim(z) <- dim(counts)
    z
}

#' @keywords internal
#' @importFrom stats ppois dpois qnorm runif
poisson_rqr <- function(counts, lambda) {
    # randomized quantile residuals for Poisson
    # The naive approach sometimes yields infinites that we can't work with
    # this splits by low and high tails and works in log space
    y <- as.vector(counts)
    lam <- as.vector(lambda)
    u <- runif(length(y))

    z <- numeric(length(y))
    bad <- is.na(lam) | is.na(y)
    low <- !bad & (y <= lam) # lower-tail side
    upp <- !bad & (y >  lam) # upper-tail side

    # lower side: p = F(y-1) + U * P(Y=y), z = qnorm(p)
    if (any(low)) {
        yl <- y[low]
        laml <- lam[low]
        lF  <- ppois(yl - 1, lambda = laml, lower.tail = TRUE, log.p = TRUE)
        ld  <- dpois(yl,     lambda = laml, log = TRUE)
        lp  <- log_add(lF, log(u[low]) + ld)
        z[low] <- qnorm(lp, log.p = TRUE)
    }

    # upper side: s = S(y) + U * P(Y=y), z = -qnorm(s)
    if (any(upp)) {
        yu <- y[upp]
        lamu <- lam[upp]
        lS  <- ppois(yu, lambda = lamu, lower.tail = FALSE, log.p = TRUE)
        ld  <- dpois(yu, lambda = lamu, log = TRUE)
        ls  <- log_add(lS, log(u[upp]) + ld)
        z[upp] <- qnorm(ls, lower.tail = FALSE, log.p = TRUE)
    }

    z[bad] <- 0
    dim(z) <- dim(counts)
    z
}

#' Apply a quantile function to standard normal draws without losing the tails
#'
#' Computes `qfun(pnorm(z), ...)` but works in log space on whichever tail is
#' closer, so that large |z| do not round to probabilities of exactly 0 or 1
#' (which would give infinite or degenerate quantiles).
#'
#' @param z Matrix of standard normal values
#' @param qfun Quantile function accepting `p`, `lower.tail` and `log.p`, e.g. `qnbinom`
#' @param params Named list of parameters for `qfun`. Each is recycled to the
#'   length of `z` (so a per-feature vector applies row-wise to a features x samples matrix).
#'
#' @return Matrix the same shape as `z`. Entries with any `NA` parameter are `NA`.
#' @keywords internal
#' @importFrom stats pnorm
quantile_from_normal <- function(z, qfun, params) {
    n <- length(z)
    zv <- as.vector(z)
    params <- lapply(params, function(p) rep_len(as.vector(p), n))
    ok <- !is.na(zv)
    for (p in params) { ok <- ok & !is.na(p) }

    out <- rep(NA_real_, n)
    low <- which(ok & zv <= 0)
    upp <- which(ok & zv > 0)
    if (length(low) > 0) {
        out[low] <- do.call(qfun, c(
            list(pnorm(zv[low], log.p = TRUE)),
            lapply(params, `[`, low),
            list(log.p = TRUE)
        ))
    }
    if (length(upp) > 0) {
        out[upp] <- do.call(qfun, c(
            list(pnorm(zv[upp], lower.tail = FALSE, log.p = TRUE)),
            lapply(params, `[`, upp),
            list(lower.tail = FALSE, log.p = TRUE)
        ))
    }
    dim(out) <- dim(z)
    out
}
