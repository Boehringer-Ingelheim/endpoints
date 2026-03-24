# util_funcs.R ---------------------------------------------------------------

# helper convention used throughout: if NULL return a, else return b
`%||%` <- function(a, b) if (!is.null(a)) a else b


# helper funcs for binary outcomes
logit <- function(p) {
  if (anyNA(p)) stop("logit(): p contains NA.")
  if (any(p <= 0 | p >= 1)) stop("logit(): p must be in (0, 1).")
  qlogis(p)
}

inv_logit <- function(x) plogis(x)

#' Solve exponential rate parameters from target observed event probabilities
#' 
#' \code{rate_from_prob()} is a helper for converting a user-specified observed
#' event probability into an exponential rate parameter under several common
#' time-to-event design settings.
#' 
#' Depending on \code{mode}, the function solves for:
#' 
#' \itemize{ \item an independent exponential censoring rate (\code{"simple"}),
#' \item an exponential event rate under pure administrative censoring
#' (\code{"admin"}), or \item an approximate censoring rate for a secondary
#' non-fatal endpoint in a semi-competing risks setting
#' (\code{"semi-competing"}). }
#' 
#' This is intended as a design-stage helper for choosing rate parameters
#' before calling \code{\link{makeData}}. Not all cases are covered. In more
#' complicated settings, we recommend iteration may be necessary.
#' 
#' 
#' @param target_prob A numeric scalar in \code{(0, 1)} giving the target
#' observed event probability.
#' 
#' The meaning depends on \code{mode}: \itemize{ \item for \code{"simple"}, the
#' desired probability that the event is observed before independent random
#' censoring; \item for \code{"admin"}, the desired probability that the event
#' occurs before the administrative follow-up limit; \item for
#' \code{"semi-competing"}, the desired approximate probability that the
#' secondary non-fatal event is observed first. }
#' @param mode Character string indicating which probability-to-rate
#' relationship to use. Must be one of: \describe{
#' \item{list("\"simple\"")}{One time-to-event endpoint with independent
#' exponential censoring.} \item{list("\"admin\"")}{One time-to-event endpoint
#' with administrative censoring only.}
#' \item{list("\"semi-competing\"")}{Approximate semi-competing risks
#' calculation for a secondary non-fatal endpoint.} }
#' @param event_rate Numeric scalar giving the exponential event rate
#' \eqn{\lambda_e}.
#' 
#' Used only when \code{mode = "simple"}.
#' @param admin_time Numeric scalar giving the administrative censoring time
#' horizon \eqn{A}.
#' 
#' Used only when \code{mode = "admin"}.
#' @param fatal_event_rate Numeric scalar giving the event rate for the primary
#' fatal event, denoted \eqn{\lambda_{e1}}.
#' 
#' Used only when \code{mode = "semi-competing"}.
#' @param fatal_censor_rate Numeric scalar giving the independent censoring
#' rate for the fatal endpoint, denoted \eqn{\lambda_{c1}}.
#' 
#' Used only when \code{mode = "semi-competing"}.
#' @param nonfatal_event_rate Numeric scalar giving the event rate for the
#' secondary non-fatal endpoint, denoted \eqn{\lambda_{e2}}.
#' 
#' Used only when \code{mode = "semi-competing"}.
#' @return A numeric scalar giving the solved rate parameter.
#' 
#' Depending on \code{mode}, this is: \itemize{ \item the censoring rate
#' (\code{"simple"}), \item the event rate (\code{"admin"}), or \item the
#' approximate non-fatal censoring rate (\code{"semi-competing"}). }
#' @section Mode = \code{"simple"}:
#' 
#' Assume a single event time \eqn{T \sim \mathrm{Exp}(\lambda_e)} and an
#' independent censoring time \eqn{C \sim \mathrm{Exp}(\lambda_c)}.
#' 
#' The probability of observing the event is
#' 
#' \deqn{ P(T < C) = \frac{\lambda_e}{\lambda_e + \lambda_c}. }
#' 
#' Given \eqn{\lambda_e} and a target probability \eqn{p}, the function solves
#' for the censoring rate \eqn{\lambda_c}:
#' 
#' \deqn{ \lambda_c = \frac{\lambda_e}{p} - \lambda_e. }
#' 
#' In this mode, the returned value is the required independent censoring rate.
#' @section Mode = \code{"admin"}:
#' 
#' Assume a single event time \eqn{T \sim \mathrm{Exp}(\lambda_e)} with fixed
#' administrative censoring at time \eqn{A > 0}.
#' 
#' The probability of observing the event by time \eqn{A} is
#' 
#' \deqn{ P(T \le A) = 1 - e^{-\lambda_e A}. }
#' 
#' Given \eqn{A} and a target probability \eqn{p}, the function solves for the
#' event rate \eqn{\lambda_e}:
#' 
#' \deqn{ \lambda_e = -\frac{\log(1 - p)}{A}. }
#' 
#' In this mode, the returned value is the required event rate.
#' @section Mode = \code{"semi-competing"}:
#' 
#' This mode uses a low-correlation / independent-clock approximation for a
#' secondary non-fatal event in a semi-competing risks setting.
#' 
#' Let: \itemize{ \item \eqn{\lambda_{e1}} be the fatal event rate, \item
#' \eqn{\lambda_{c1}} be the censoring rate for the fatal endpoint, \item
#' \eqn{\lambda_{e2}} be the non-fatal event rate, \item \eqn{\lambda_{c2}} be
#' the censoring rate for the non-fatal endpoint. }
#' 
#' Under the approximation used in the package vignette, the probability that
#' the secondary non-fatal event is observed first is
#' 
#' \deqn{ p \approx \frac{\lambda_{e2}} {\lambda_{e1} + \lambda_{c1} +
#' \lambda_{e2} + \lambda_{c2}}. }
#' 
#' Given \eqn{p}, \eqn{\lambda_{e1}}, \eqn{\lambda_{c1}}, and
#' \eqn{\lambda_{e2}}, the function solves for \eqn{\lambda_{c2}}:
#' 
#' \deqn{ \lambda_{c2} = \frac{\lambda_{e2}}{p} - \left(\lambda_{e1} +
#' \lambda_{c1} + \lambda_{e2}\right). }
#' 
#' In this mode, the returned value is the approximate non-fatal censoring
#' rate.
#' 
#' Because this is an approximation, the realized observed non-fatal event
#' probability in \code{\link{makeData}} may differ once dependence, censoring,
#' and fatal-event logic are fully applied.
#' @seealso \code{\link{makeData}} for simulation of trial datasets using the
#' resulting rates.
#' @examples
#' 
#' # ------------------------------------------------------------
#' # 1) Independent exponential censoring only
#' # Solve for censoring rate given event rate and target event probability
#' # ------------------------------------------------------------
#' rate_from_prob(
#'   target_prob = 0.90,
#'   mode = "simple",
#'   event_rate = 1 / 24
#' )
#' # approximately 1/216
#' 
#' 
#' # ------------------------------------------------------------
#' # 2) Administrative censoring only
#' # Solve for event rate given follow-up time and target event probability
#' # ------------------------------------------------------------
#' rate_from_prob(
#'   target_prob = 0.20,
#'   mode = "admin",
#'   admin_time = 4
#' )
#' # approximately 0.0558
#' 
#' 
#' # ------------------------------------------------------------
#' # 3) Semi-competing risks approximation
#' # Solve for the non-fatal censoring rate
#' # ------------------------------------------------------------
#' rate_from_prob(
#'   target_prob = 0.20,
#'   mode = "semi-competing",
#'   fatal_event_rate = 1 / 50,
#'   fatal_censor_rate = 1 / 16.667,
#'   nonfatal_event_rate = 1 / 35
#' )
#' 
#' @export
rate_from_prob <- function(target_prob,
                           mode = c("simple", "admin", "semi-competing"),
                           event_rate = NULL,
                           admin_time = NULL,
                           fatal_event_rate = NULL,
                           fatal_censor_rate = NULL,
                           nonfatal_event_rate = NULL) {

  mode <- match.arg(mode)

  if (!is.numeric(target_prob) || length(target_prob) != 1L ||
      is.na(target_prob) || target_prob <= 0 || target_prob >= 1) {
    stop("`target_prob` must be a single number in (0, 1).")
  }

  # ---- simple: one endpoint + independent exponential censoring ----------
  if (mode == "simple") {
    if (!is.numeric(event_rate) || length(event_rate) != 1L ||
        is.na(event_rate) || event_rate <= 0) {
      stop("For mode = 'simple', `event_rate` must be a single positive number.")
    }

    return(event_rate / target_prob - event_rate)
  }

  # ---- admin: administrative censoring only ------------------------------
  if (mode == "admin") {
    if (!is.numeric(admin_time) || length(admin_time) != 1L ||
        is.na(admin_time) || admin_time <= 0) {
      stop("For mode = 'admin', `admin_time` must be a single positive number.")
    }

    return(-log(1 - target_prob) / admin_time)
  }

  # ----  semi-competing risks approximation ---------------------------
  # p ~= lambda_e2 / (lambda_e1 + lambda_c1 + lambda_e2 + lambda_c2)
  # solve for lambda_c2
  if (mode == "semi-competing") {
    vals <- c(fatal_event_rate, fatal_censor_rate, nonfatal_event_rate)
    if (!is.numeric(vals) || anyNA(vals) || any(vals < 0)) {
      stop(
        "For mode = 'src', `fatal_event_rate`, `fatal_censor_rate`, and ",
        "`nonfatal_event_rate` must be numeric and >= 0."
      )
    }
    if (length(fatal_event_rate) != 1L ||
        length(fatal_censor_rate) != 1L ||
        length(nonfatal_event_rate) != 1L) {
      stop(
        "For mode = 'src', `fatal_event_rate`, `fatal_censor_rate`, and ",
        "`nonfatal_event_rate` must each be length 1."
      )
    }

    out <- nonfatal_event_rate / target_prob -
      (fatal_event_rate + fatal_censor_rate + nonfatal_event_rate)

    if (out < 0) {
      warning(
        "Computed censoring rate is negative; `target_prob` may be infeasible ",
        "under the semi-competing risks approximation."
      )
    }

    return(out)
  }
}

#' Construct a correlation matrix from endpoint-index triplets
#' 
#' \code{corr_make()} constructs a correlation matrix from a set of \eqn{(i, j,
#' \rho)} triplets, where \eqn{i} and \eqn{j} are endpoint indices and
#' \eqn{\rho} is their desired pairwise correlation.
#' 
#' This is a convenience function for specifying sparse correlation structures
#' when using \code{\link{makeData}}. Any unspecified off-diagonal entries are
#' left at 0, and diagonal entries are set to 1.
#' 
#' 
#' @param num_endpoints A single positive integer giving the number of
#' endpoints (i.e., the dimension of the correlation matrix to construct).
#' @param values Optional object coercible to a numeric matrix with 3 columns.
#' Each row must have the form \code{c(i, j, rho)}, where: \describe{
#' \item{list("i")}{First endpoint index.} \item{list("j")}{Second endpoint
#' index.} \item{list("rho")}{Desired correlation between endpoints \code{i}
#' and \code{j}.} }
#' 
#' \code{values} may be a matrix, data frame, or vector coercible to 3 columns.
#' 
#' The endpoint indices in values correspond to the order in which endpoints
#' are supplied in endpoint_details; for example, index 1 refers to the first
#' endpoint in the list, index 2 to the second, and so on.
#' 
#' If \code{NULL}, the identity matrix of size \code{num_endpoints} is
#' returned. Unspecified off-diagonal entries default to 0.
#' @return A symmetric numeric correlation matrix of dimension
#' \code{num_endpoints x num_endpoints}, with ones on the diagonal.
#' @section How entries are interpreted:
#' 
#' For each row in \code{values}, the indices \code{i} and \code{j} are mapped
#' to the upper triangle of the matrix using: \deqn{ i^\star = \min(i, j),
#' \qquad j^\star = \max(i, j), } so the order of the indices does not matter.
#' The specified correlation is then assigned symmetrically: \deqn{ R[i^\star,
#' j^\star] = R[j^\star, i^\star] = \rho. } Finally, the diagonal is reset to
#' exactly 1, so any user-supplied diagonal entries are overwritten.
#' @section Use with \code{makeData()}:
#' 
#' The resulting matrix can be passed directly to \code{\link{makeData}} as
#' \code{correlation_matrix}.
#' 
#' When \code{target_correlation = FALSE}, the resulting matrix is interpreted
#' as the latent Gaussian correlation matrix.
#' 
#' When \code{target_correlation = TRUE}, the resulting matrix is treated as
#' the desired observed Pearson correlation matrix, and the simulation routine
#' attempts to calibrate a latent Gaussian matrix to approximately match it
#' after marginal transformation.
#' @seealso \code{\link{makeData}} for the main simulation function.
#' 
#' \code{\link{endpoint_details}} for endpoint specification details.
#' 
#' \code{\link{calibration_control}} for correlation calibration settings.
#' @keywords utilities
#' @examples
#' 
#' ## Identity correlation matrix
#' corr_make(num_endpoints = 3)
#' 
#' ## Three endpoints with pairwise specification
#' R3 <- corr_make(
#'   num_endpoints = 3,
#'   values = rbind(
#'     c(1, 2, 0.20),
#'     c(1, 3, 0.10),
#'     c(2, 3, 0.15)
#'   )
#' )
#' R3
#' 
#' ## AR(1)-style structure across five repeated measures
#' rho <- 0.5
#' R5 <- corr_make(
#'   num_endpoints = 5,
#'   values = rbind(
#'     c(1, 2, rho),   c(1, 3, rho^2), c(1, 4, rho^3), c(1, 5, rho^4),
#'     c(2, 3, rho),   c(2, 4, rho^2), c(2, 5, rho^3),
#'     c(3, 4, rho),   c(3, 5, rho^2),
#'     c(4, 5, rho)
#'   )
#' )
#' R5
#' 
#' ## Order of indices does not matter
#' corr_make(
#'   num_endpoints = 3,
#'   values = rbind(
#'     c(2, 1, 0.30),
#'     c(3, 1, 0.10)
#'   )
#' )
#'
#' @export
corr_make <- function(num_endpoints, values = NULL) {
  # ---- validate dimension --------------------------------------------------
  if (!is.numeric(num_endpoints) || length(num_endpoints) != 1L || is.na(num_endpoints)) {
    stop("`num_endpoints` must be a single positive integer.")
  }
  num_endpoints <- as.integer(num_endpoints)
  if (num_endpoints < 1L) {
    stop("`num_endpoints` must be >= 1.")
  }

  # start from identity correlation matrix
  R <- diag(num_endpoints)

  # nothing else to fill
  if (is.null(values)) {
    return(R)
  }

  # ---- coerce and validate entries ----------------------------------------
  vals <- matrix(as.matrix(values), ncol = 3)

  if (nrow(vals) == 0L) {
    return(R)
  }

  ij  <- vals[, 1:2, drop = FALSE]
  rho <- vals[, 3]

  if (anyNA(vals)) {
    stop("`values` must not contain NA.")
  }
  if (any(ij <= 0)) {
    stop("Indices in `values[, 1:2]` must be positive.")
  }
  if (any(ij > num_endpoints)) {
    stop("Indices in `values[, 1:2]` exceed `num_endpoints`.")
  }
  if (any(rho < -1 | rho > 1)) {
    stop("All correlations in `values[, 3]` must be in [-1, 1].")
  }

  # ---- map to upper triangle, assign, symmetrize ---------------------------
  i <- pmin.int(ij[, 1], ij[, 2])
  j <- pmax.int(ij[, 1], ij[, 2])

  R[cbind(i, j)] <- rho
  R[cbind(j, i)] <- rho   # symmetric fill in one shot

  # enforce correlation diagonal exactly
  diag(R) <- 1

  R
}


normalize_endpoint_type <- function(x) {
  # Aliases for use in makeData(), endpoint_details
  x <- tolower(x)
  if (x %in% c("normal", "continuous", "gaussian")) return("continuous")
  if (x %in% c("binary", "bin")) return("binary")
  if (x %in% c("count", "nb", "negbinom", "negative binomial", "zinb")) return("count")
  if (x %in% c("tte", "time-to-event", "time_to_event", "survival")) return("time-to-event")
  stop("Unknown endpoint_type: ", x)
}


# ensure correlation matrix is PD
make_pd <- function(R, eps = 1e-10) {
  out <- tryCatch(chol(R), error = function(e) NULL)
  if (!is.null(out)) return(list(R = R, L = out))

  ee <- eigen(R, symmetric = TRUE)
  vals2 <- pmax(ee$values, eps)
  R2 <- ee$vectors %*% diag(vals2, nrow = length(vals2)) %*% t(ee$vectors)

  d <- sqrt(diag(R2))
  R2 <- sweep(sweep(R2, 1, d, "/"), 2, d, "/")

  L2 <- chol(R2)
  list(R = R2, L = L2)
}


# Project to nearest positive-definite correlation matrix.
# Prefer Matrix::nearPD if available; otherwise fall back to eigen-based fix.
project_to_pd_cor <- function(R, conv_norm_type = "F") {
  if (requireNamespace("Matrix", quietly = TRUE)) {
    pmt <- Matrix::nearPD(
      R, corr = TRUE, keepDiag = TRUE,
      conv.norm.type = conv_norm_type, trace = FALSE
    )
    return(as.matrix(pmt$mat))
  }
  make_pd(R)$R
}

# -----------------------------------------------------------------------------
# Gaussian copula
# -----------------------------------------------------------------------------

# Generate an n x p matrix of dependent Uniform(0,1) variables with Gaussian copula
# dependence using Cholesky factor L of the latent correlation matrix.
draw_gaussian_copula_u <- function(n, L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- ncol(L)
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% L # similar to MASS::mvrnorm()
  U <- pnorm(Z)
  eps <- 1e-12
  pmin(pmax(U, eps), 1 - eps)
}

# -----------------------------------------------------------------------------
# Marginal quantile helpers
# -----------------------------------------------------------------------------

#  zero-inflated negative binomial quantile function:
# with prob p0 -> 0; else NB(mu, size)
# avoids dependence on external pkg
qzinb_mixture <- function(u, mu, size, p0) {
  out <- integer(length(u))
  is0 <- (u < p0)
  out[is0] <- 0L
  u2 <- (u[!is0] - p0) / (1 - p0)
  out[!is0] <- qnbinom(u2, size = size, mu = mu)
  out
}


# -----------------------------------------------------------------------------
# Correlation calibration (empirical)
# -----------------------------------------------------------------------------

# Pairwise latent-correlation calibration:
# Find rho_Z such that cor(Q1(Phi(Z1)), Q2(Phi(Z2))) ~= target_r
calibrate_latent_rho_pair <- function(
    target_r,
    qfun1,
    qfun2,
    n_mc    = 10000,
    seed    = NULL,
    tol     = 0.001,
    maxit   = 100,
    rho_cap = 0.999
) {
  # ---- checks --------------------------------------------------------------
  if (!is.numeric(target_r) || length(target_r) != 1L || is.na(target_r)) {
    stop("`target_r` must be a single numeric value.")
  }
  if (target_r < -1 || target_r > 1) {
    stop("`target_r` must be in [-1, 1].")
  }
  if (!is.function(qfun1) || !is.function(qfun2)) {
    stop("`qfun1` and `qfun2` must be functions.")
  }
  if (!is.numeric(n_mc) || length(n_mc) != 1L || is.na(n_mc) || n_mc <= 0) {
    stop("`n_mc` must be a positive number.")
  }
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("`tol` must be > 0.")
  }
  if (!is.numeric(maxit) || length(maxit) != 1L || is.na(maxit) || maxit <= 0) {
    stop("`maxit` must be > 0.")
  }
  if (!is.numeric(rho_cap) || length(rho_cap) != 1L || is.na(rho_cap) ||
      rho_cap <= 0 || rho_cap >= 1) {
    stop("`rho_cap` must be in (0, 1).")
  }

  target_r <- max(min(target_r, rho_cap), -rho_cap)

  if (!is.null(seed)) set.seed(seed)
  z1_base <- rnorm(n_mc)
  z2_base <- rnorm(n_mc)

  observed_corr_minus_target <- function(rho) {
    rho <- max(min(rho, rho_cap), -rho_cap)

    z1 <- z1_base
    z2 <- rho * z1_base + sqrt(pmax(0, 1 - rho^2)) * z2_base

    u1 <- pnorm(z1)
    u2 <- pnorm(z2)

    x1 <- qfun1(u1)
    x2 <- qfun2(u2)

    if (sd(x1) == 0 || sd(x2) == 0) return(NA_real_)
    cor(x1, x2) - target_r
  }

  lo <- -rho_cap
  hi <-  rho_cap
  f_lo <- observed_corr_minus_target(lo)
  f_hi <- observed_corr_minus_target(hi)

  if (is.na(f_lo) || is.na(f_hi)) {
    return(list(root = 0, note = "degenerate margin"))
  }

  # infeasible target after transform -> choose closest endpoint
  if (f_lo * f_hi > 0) {
    root_pick <- if (abs(f_lo) < abs(f_hi)) lo else hi
    return(list(root = root_pick, note = "target infeasible (picked closest endpoint)"))
  }

  stats::uniroot(
    f = observed_corr_minus_target,
    interval = c(lo, hi),
    tol = tol,
    maxiter = maxit
  )
}


# Matrix-level latent correlation calibration for Gaussian copula simulation.
# Returns latent correlation matrix to approximately match target observed Pearson corr
# after applying endpoint-specific quantile maps.
calibrate_latent_cor_matrix <- function(
    target_cor,
    qfuns,
    ensure_pd = TRUE,
    conv_norm_type = "F",
    return_diagnostics = FALSE,
    n_mc    = 10000,
    seed    = NULL,
    tol     = 0.001,
    maxit   = 100,
    rho_cap = 0.999
) {
  # ---- checks --------------------------------------------------------------
  if (!is.matrix(target_cor) || nrow(target_cor) != ncol(target_cor)) {
    stop("`target_cor` must be a square matrix.")
  }
  if (max(abs(target_cor - t(target_cor))) > 1e-10) {
    stop("`target_cor` must be symmetric.")
  }
  if (max(abs(diag(target_cor) - 1)) > 1e-10) {
    stop("`target_cor` must have 1s on the diagonal.")
  }
  if (any(target_cor < -1 | target_cor > 1, na.rm = TRUE)) {
    stop("All entries of `target_cor` must be in [-1, 1].")
  }

  p <- ncol(target_cor)

  if (!is.list(qfuns) || length(qfuns) != p) {
    stop("`qfuns` must be a list of length ncol(target_cor).")
  }
  if (!all(vapply(qfuns, is.function, logical(1)))) {
    stop("All elements of `qfuns` must be functions.")
  }

  if (!is.logical(ensure_pd) || length(ensure_pd) != 1L) {
    stop("`ensure_pd` must be TRUE/FALSE.")
  }
  if (!is.logical(return_diagnostics) || length(return_diagnostics) != 1L) {
    stop("`return_diagnostics` must be TRUE/FALSE.")
  }

  latent_cor <- diag(p)
  conv_diag <- vector("list", p)

  for (r in seq_len(p - 1L)) {
    conv_diag[[r]] <- vector("list", p)

    for (c in (r + 1L):p) {
      seed_rc <- if (is.null(seed)) NULL else (seed + 1000L * r + c)

      fit_rc <- calibrate_latent_rho_pair(
        target_r = target_cor[r, c],
        qfun1    = qfuns[[r]],
        qfun2    = qfuns[[c]],
        n_mc     = n_mc,
        seed     = seed_rc,
        tol      = tol,
        maxit    = maxit,
        rho_cap  = rho_cap
      )

      latent_cor[r, c] <- fit_rc$root
      conv_diag[[r]][[c]] <- fit_rc
    }
  }

  latent_cor[lower.tri(latent_cor)] <- t(latent_cor)[lower.tri(latent_cor)]

  if (isTRUE(ensure_pd)) {
    latent_cor <- project_to_pd_cor(latent_cor, conv_norm_type = conv_norm_type)
  }

  if (isTRUE(return_diagnostics)) {
    return(list(cor_mat = latent_cor, convergence = conv_diag))
  }

  latent_cor
}

# -----------------------------------------------------------------------------
# S3 constructor + methods
# -----------------------------------------------------------------------------

new_makeDataSim <- function(data, meta) {
  stopifnot(is.data.frame(data))
  stopifnot(is.list(meta))
  structure(
    list(data = data, meta = meta),
    class = "makeDataSim"
  )
}

#' @exportS3Method
as.data.frame.makeDataSim <- function(x, ...) {
  x$data
}

#' @exportS3Method
print.makeDataSim <- function(x, ...) {
  dat  <- x$data
  meta <- x$meta

  has_trt <- "trt" %in% names(dat)
  K <- meta$n_arms %||% if (has_trt) length(unique(dat$trt)) else 1L

  cat("<makeDataSim>\n")
  cat("  n =", nrow(dat), "\n", sep = "")
  cat("  n_arms =", K, "\n", sep = "")
  cat("  endpoints =", length(meta$endpoint_details %||% list()), "\n", sep = "")
  cat("  target_correlation =", isTRUE(meta$target_correlation), "\n", sep = "")
  cat("  single_endpoint_mode =", isTRUE(meta$single_endpoint_mode), "\n", sep = "")
  cat("  control_only =", isTRUE(meta$control_only), "\n", sep = "")

  if (has_trt) {
    arm_tab <- table(dat$trt)
    cat("  n_by_arm = ")
    cat(paste0(names(arm_tab), ":", as.integer(arm_tab), collapse = ", "))
    cat("\n")
  } else {
    cat("  n_by_arm = 0:", nrow(dat), "\n", sep = "")
  }

  if (!is.null(meta$endpoint_types)) {
    cat("  endpoint_types = ", paste(meta$endpoint_types, collapse = ", "), "\n", sep = "")
  }

  show_rows <- 6L
  show_cols <- 10L

  disp <- dat
  n <- nrow(disp); p <- ncol(disp)

  if (p > show_cols) {
    keep <- c(seq_len(show_cols - 1L), p)
    disp <- disp[, keep, drop = FALSE]
    names(disp)[show_cols] <- paste0("\u2026", names(disp)[show_cols])
  }

  if (n > show_rows) {
    disp_head <- disp[seq_len(show_rows - 1L), , drop = FALSE]
    ell_row <- as.list(rep("\u22ee", ncol(disp)))
    names(ell_row) <- names(disp)
    disp <- rbind(disp_head, as.data.frame(ell_row, stringsAsFactors = FALSE))
    rownames(disp) <- c(rownames(disp_head), "\u22ee")
  }

  cat("\n  data (head):\n")
  print(disp)

  cat("\n")
  invisible(x)
}


#' Summarize simulated data from a \code{makeDataSim} object
#' 
#' \code{summary.makeDataSim()} computes arm-specific diagnostic summaries for
#' a \code{"makeDataSim"} object returned by \code{\link{makeData}}.
#' 
#' The summary includes: \itemize{ \item the target correlation matrix supplied
#' to \code{makeData()}, \item empirical endpoint correlations within each
#' treatment arm, \item endpoint-specific marginal summaries for continuous,
#' binary, count, and time-to-event outcomes. }
#' 
#' These summaries are intended as a quick validation tool for checking that
#' the simulated dataset is broadly consistent with the requested
#' data-generating parameters.
#' 
#' 
#' @param object An object of class \code{"makeDataSim"}, typically created by
#' \code{\link{makeData}}.
#' @param \dots Currently unused. Included for S3 method compatibility.
#' @return A list of class \code{"summary.makeDataSim"} with components:
#' \describe{ \item{list("target_correlation")}{The target correlation matrix
#' supplied to \code{makeData()}, or \code{NULL} in single-endpoint mode.}
#' \item{list("estimated_correlation_by_arm")}{A named list of empirical
#' arm-specific correlation matrices.} \item{list("continuous")}{A data frame
#' of continuous-endpoint summaries, or \code{NULL} if no continuous endpoints
#' were simulated.} \item{list("binary")}{A data frame of binary-endpoint
#' summaries, or \code{NULL} if no binary endpoints were simulated.}
#' \item{list("count")}{A data frame of count-endpoint summaries, or
#' \code{NULL} if no count endpoints were simulated.} \item{list("tte")}{A data
#' frame of time-to-event summaries, or \code{NULL} if no TTE endpoints were
#' simulated.} \item{list("n_arms")}{The total number of treatment arms
#' represented in the simulated dataset.} }
#' @section General behavior:
#' 
#' The summary method extracts the simulated dataset and metadata stored in the
#' \code{"makeDataSim"} object, then computes: \enumerate{ \item the requested
#' target correlation matrix (as stored in
#' \code{object$meta$correlation_matrix}), \item the observed Pearson
#' correlation matrix among endpoint columns within each treatment arm, \item
#' endpoint-type-specific summaries based on simple fitted models or direct
#' descriptive estimators. }
#' @section Arm-specific empirical correlation:
#' 
#' For each study arm, the method computes the Pearson correlation matrix of
#' the simulated endpoint columns listed in \code{object$meta$endpoint_names}.
#' These are returned as a named list in \code{$estimated_correlation_by_arm},
#' with elements named \code{"arm_0"}, \code{"arm_1"}, and so on.
#' 
#' Correlations are computed using \code{cor()} and rounded to 3 decimal
#' places. %If some endpoint combinations are degenerate or discrete, warnings
#' from \code{cor()} are suppressed.
#' @section Continuous endpoints:
#' 
#' For each continuous endpoint: \itemize{ \item a linear model is fit
#' (\code{lm(y ~ trt)} in multi-arm settings, otherwise \code{lm(y ~ 1)}),
#' \item the control-group intercept is used as the estimated baseline mean,
#' \item treatment coefficients are reported as estimated mean shifts, \item
#' arm-specific residual SDs are computed from the residuals of the fitted
#' model. }
#' 
#' The returned table includes: \describe{ \item{list("endpoint")}{Endpoint
#' name (e.g. \code{Cont_1}).} \item{list("arm")}{Arm index.}
#' \item{list("input_baseline_mean")}{Encoded control-group mean.}
#' \item{list("input_sd")}{Encoded SD for that arm.}
#' \item{list("input_trt_effect")}{Encoded treatment effect for that arm.}
#' \item{list("est_baseline_mean")}{Estimated control-group mean from the
#' fitted model.} \item{list("est_trt_effect")}{Estimated mean shift for that
#' arm vs control.} \item{list("est_resid_sd")}{Residual SD within that arm.} }
#' @section Binary endpoints:
#' 
#' For each binary endpoint: \itemize{ \item a logistic regression model is fit
#' (\code{glm(..., family = binomial())}), \item the control-group intercept is
#' interpreted on the logit scale, \item treatment coefficients are reported as
#' estimated log-odds ratios, \item observed arm-specific event probabilities
#' are computed directly as sample means. }
#' 
#' The returned table includes: \describe{ \item{list("endpoint")}{Endpoint
#' name (e.g. \code{Bin_1}).} \item{list("arm")}{Arm index.}
#' \item{list("input_baseline_prob")}{Encoded control-group probability.}
#' \item{list("input_trt_logOR")}{Encoded treatment effect on the log-odds
#' scale.} \item{list("input_trt_prob")}{Encoded arm-specific probability.}
#' \item{list("est_baseline_prob")}{Estimated control-group probability from
#' the fitted model.} \item{list("est_trt_logOR")}{Estimated treatment log-odds
#' ratio for that arm vs control.} \item{list("est_prob")}{Observed
#' arm-specific event proportion.} }
#' @section Count endpoints:
#' 
#' For each count endpoint: \itemize{ \item a negative binomial model is fit
#' using \code{MASS::glm.nb()}, \item the control-group intercept is
#' exponentiated to obtain the estimated baseline mean, \item treatment
#' coefficients are reported as estimated log rate-ratios, \item the fitted
#' dispersion parameter is reported, \item observed arm-specific means and
#' structural zero proportions are computed directly. }
#' 
#' The returned table includes: \describe{ \item{list("endpoint")}{Endpoint
#' name (e.g. \code{Int_1}).} \item{list("arm")}{Arm index.}
#' \item{list("input_baseline_mean")}{Specified control-group mean.}
#' \item{list("input_trt_logRR")}{Specified treatment effect on the log
#' rate-ratio scale.} \item{list("input_trt_mean")}{Specified arm-specific mean
#' count.} \item{list("input_size")}{Specified negative-binomial size
#' parameter.} \item{list("input_p_zero")}{Specified structural zero
#' probability.} \item{list("est_baseline_mean")}{Estimated control-group mean
#' from the fitted model.} \item{list("est_trt_logRR")}{Estimated treatment log
#' rate-ratio for that arm vs control.} \item{list("est_size")}{Estimated
#' negative-binomial size parameter.} \item{list("obs_mean")}{Observed
#' arm-specific mean count.} \item{list("obs_p0")}{Observed proportion of zeros
#' in that arm.} }
#' @section Time-to-event endpoints:
#' 
#' For each time-to-event endpoint: \itemize{ \item a Cox proportional hazards
#' model is fit using \code{survival::coxph()}, \item treatment coefficients
#' are reported as estimated log hazard-ratios, \item observed arm-specific
#' event proportions are reported, \item an arm-specific exponential maximum
#' likelihood estimate of the event rate is computed as \eqn{\hat\lambda =
#' \sum_i \delta_i / \sum_i t_i}, where \eqn{t_i} is observed follow-up time
#' and \eqn{\delta_i} is the event indicator. }
#' 
#' The returned table includes: \describe{ \item{list("endpoint")}{Endpoint
#' name (e.g. \code{TTE_1}).} \item{list("arm")}{Arm index.}
#' \item{list("censor_col")}{Corresponding event-indicator column name (e.g.
#' \code{Status_1}).} \item{list("input_baseline_rate")}{Requested
#' control-group exponential event rate.}
#' \item{list("input_trt_logHR")}{Requested treatment effect on the log
#' hazard-ratio scale.} \item{list("est_trt_logHR")}{Estimated log hazard-ratio
#' from the Cox model.} \item{list("input_trt_HR")}{Requested treatment effect
#' on the hazard-ratio scale.} \item{list("est_trt_HR")}{Estimated hazard-ratio
#' from the Cox model.} \item{list("obs_event_rate")}{Observed event proportion
#' in that arm.} \item{list("exp_rate")}{Arm-specific exponential MLE
#' event-rate estimate.} }
#' @section Interpretation:
#' 
#' Because simulations are finite-sample and may involve nonlinear marginal
#' transformations, censoring, administrative censoring, and fatal/non-fatal
#' event logic, the estimated summaries will generally not match the requested
#' inputs exactly. The summary output is therefore best viewed as a diagnostic
#' check rather than an exact recovery target.
#' @seealso \code{\link{makeData}} for the main simulation routine.
#' 
#' \code{\link[=print.summary.makeDataSim]{print.summary.makeDataSim}} for the
#' print method.
#' 
#' \code{\link{plot.makeDataSim}} for graphical diagnostics.
#' 
#' \code{\link{corr_make}} for building correlation matrices.
#' @keywords methods
#' @examples
#' 
#' ## Continuous + binary example
#' ep1 <- list(
#'   endpoint_type = "continuous",
#'   baseline_mean = 10,
#'   sd            = 2,
#'   trt_effect    = -1
#' )
#' 
#' ep2 <- list(
#'   endpoint_type = "binary",
#'   baseline_prob = 0.30,
#'   trt_prob      = 0.45
#' )
#' 
#' R2 <- corr_make(
#'   num_endpoints = 2,
#'   values = rbind(c(1, 2, 0.2))
#' )
#' 
#' sim_obj <- makeData(
#'   correlation_matrix    = R2,
#'   sample_size_per_group = 500,
#'   SEED                  = 1,
#'   endpoint_details      = list(ep1, ep2)
#' )
#' 
#' ss <- summary(sim_obj)
#' 
#' ## Display top-level structure
#' ss$n_arms
#' ss$target_correlation
#' ss$estimated_correlation_by_arm
#' 
#' ## Endpoint-specific summaries
#' ss$continuous
#' ss$binary
#' 
#' ## Time-to-event example
#' ep_tte <- list(
#'   endpoint_type  = "tte",
#'   baseline_rate  = 1 / 24,
#'   trt_effect     = log(0.8),
#'   censoring_rate = 1 / 216,
#'   fatal_event    = TRUE
#' )
#' 
#' sim_tte <- makeData(
#'   correlation_matrix    = NULL,
#'   sample_size_per_group = 1000,
#'   SEED                  = 2,
#'   endpoint_details      = list(ep_tte)
#' )
#' 
#' summary(sim_tte)$tte
#' 
#' @exportS3Method
summary.makeDataSim <- function(object, ...) {
  dat  <- object$data
  meta <- object$meta

  endpoint_names   <- meta$endpoint_names
  endpoint_types   <- meta$endpoint_types
  endpoint_details <- meta$endpoint_details
  cor_target       <- meta$correlation_matrix

  has_trt <- "trt" %in% names(dat)
  K <- meta$n_arms %||% if (has_trt) length(unique(dat$trt)) else 1L

  trtF <- if (has_trt) {
    factor(dat$trt, levels = 0:(K - 1L))
  } else {
    factor(rep(0, nrow(dat)), levels = 0)
  }

  cor_by_arm <- setNames(vector("list", K), paste0("arm_", 0:(K - 1L)))
  for (a in 0:(K - 1L)) {
    Xa <- if (has_trt) dat[dat$trt == a, endpoint_names, drop = FALSE] else dat[, endpoint_names, drop = FALSE]
    cor_by_arm[[paste0("arm_", a)]] <- round(suppressWarnings(cor(Xa)), 3)
  }

  idx_cont <- which(endpoint_types == "continuous")
  idx_bin  <- which(endpoint_types == "binary")
  idx_cnt  <- which(endpoint_types == "count")
  idx_tte  <- which(endpoint_types == "time-to-event")

  expand_input_effect <- function(eff, K) {
    if (is.null(eff)) return(rep(0, K))
    if (length(eff) == 1L) return(c(0, rep(eff, K - 1L)))
    if (length(eff) == (K - 1L)) return(c(0, eff))
    stop("Internal: invalid trt_effect length; expected NULL, 1, or K-1.")
  }

  expand_active <- function(x, K) {
    if (is.null(x)) return(NULL)
    if (length(x) == 1L) return(rep(x, K - 1L))
    if (length(x) == (K - 1L)) return(as.numeric(x))
    stop("Internal: invalid active-arm length; expected 1 or K-1.")
  }

  arm_subset <- function(a) {
    if (has_trt) dat$trt == a else rep(TRUE, nrow(dat))
  }

  cont_tbl <- NULL
  if (length(idx_cont) > 0) {
    cont_tbl <- do.call(rbind, lapply(idx_cont, function(j) {
      col  <- endpoint_names[j]
      spec <- endpoint_details[[j]]

      fit <- if (has_trt && K > 1L) lm(dat[[col]] ~ trtF) else lm(dat[[col]] ~ 1)
      b0 <- unname(coef(fit)[1])

      s0 <- unname(spec$baseline_mean)
      sd_in <- spec$sd

      if (length(sd_in) == 1L) {
        sd_in_vec <- rep(sd_in, K)
      } else if (length(sd_in) == K) {
        sd_in_vec <- as.numeric(sd_in)
      } else {
        sd_in_vec <- rep(NA_real_, K)
      }

      est_eff <- rep(0, K)
      names(est_eff) <- 0:(K - 1L)

      if (has_trt && K > 1L) {
        cf <- coef(fit)
        for (a in 1:(K - 1L)) {
          nm <- paste0("trtF", a)
          if (nm %in% names(cf)) est_eff[a + 1L] <- unname(cf[nm])
        }
      }

      inp_eff <- expand_input_effect(spec$trt_effect %||% NULL, K)

      est_sd_by_arm <- vapply(0:(K - 1L), function(a) {
        idx <- arm_subset(a)
        if (sum(idx) <= 1) return(NA_real_)
        stats::sd(residuals(fit)[idx])
      }, numeric(1))

      do.call(rbind, lapply(0:(K - 1L), function(a) {
        data.frame(
          endpoint = col,
          arm = a,
          input_baseline_mean = s0,
          input_sd            = sd_in_vec[a + 1L],
          input_trt_effect    = inp_eff[a + 1L],
          est_baseline_mean   = b0,
          est_trt_effect      = est_eff[a + 1L],
          est_resid_sd        = est_sd_by_arm[a + 1L],
          stringsAsFactors = FALSE
        )
      }))
    }))
    rownames(cont_tbl) <- NULL
  }

  bin_tbl <- NULL
  if (length(idx_bin) > 0) {
    bin_tbl <- do.call(rbind, lapply(idx_bin, function(j) {
      col  <- endpoint_names[j]
      spec <- endpoint_details[[j]]

      fit <- if (has_trt && K > 1L) {
        glm(dat[[col]] ~ trtF, family = binomial())
      } else {
        glm(dat[[col]] ~ 1, family = binomial())
      }

      cf <- coef(fit)
      b0 <- unname(cf[1])

      est_eff <- rep(0, K)
      names(est_eff) <- 0:(K - 1L)
      if (has_trt && K > 1L) {
        for (a in 1:(K - 1L)) {
          nm <- paste0("trtF", a)
          if (nm %in% names(cf)) est_eff[a + 1L] <- unname(cf[nm])
        }
      }

      inp_eff <- expand_input_effect(spec$trt_effect %||% NULL, K)
      inp_prob_active <- expand_active(spec$trt_prob %||% NULL, K)
      inp_prob_vec <- rep(spec$baseline_prob, K)
      if (!is.null(inp_prob_active) && K > 1L) {
        inp_prob_vec[2:K] <- inp_prob_active
        inp_eff <- c(0, logit(inp_prob_active) - logit(spec$baseline_prob))
      }

      est_prob_by_arm <- vapply(0:(K - 1L), function(a) {
        idx <- arm_subset(a)
        mean(dat[[col]][idx])
      }, numeric(1))

      do.call(rbind, lapply(0:(K - 1L), function(a) {
        data.frame(
          endpoint = col,
          arm = a,
          input_baseline_prob = spec$baseline_prob,
          input_trt_logOR     = inp_eff[a + 1L],
          input_trt_prob      = inp_prob_vec[a + 1L],
          est_baseline_prob   = inv_logit(b0),
          est_trt_logOR       = est_eff[a + 1L],
          est_prob            = est_prob_by_arm[a + 1L],
          stringsAsFactors = FALSE
        )
      }))
    }))
    rownames(bin_tbl) <- NULL
  }

  cnt_tbl <- NULL
  if (length(idx_cnt) > 0) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("MASS is required for count endpoint summary via MASS::glm.nb().")
    }

    cnt_tbl <- do.call(rbind, lapply(idx_cnt, function(j) {
      col  <- endpoint_names[j]
      spec <- endpoint_details[[j]]

      fit <- if (has_trt && K > 1L) {
        MASS::glm.nb(dat[[col]] ~ trtF)
      } else {
        MASS::glm.nb(dat[[col]] ~ 1)
      }

      cf <- coef(fit)
      b0 <- unname(cf[1])

      est_eff <- rep(0, K)
      names(est_eff) <- 0:(K - 1L)
      if (has_trt && K > 1L) {
        for (a in 1:(K - 1L)) {
          nm <- paste0("trtF", a)
          if (nm %in% names(cf)) est_eff[a + 1L] <- unname(cf[nm])
        }
      }

      inp_eff <- expand_input_effect(spec$trt_effect %||% NULL, K)
      inp_mean_active <- expand_active(spec$trt_count %||% NULL, K)
      inp_mean_vec <- rep(spec$baseline_mean, K)
      if (!is.null(inp_mean_active) && K > 1L) {
        inp_mean_vec[2:K] <- inp_mean_active
        inp_eff <- c(0, log(inp_mean_active / spec$baseline_mean))
      }

      obs_mean_by_arm <- vapply(0:(K - 1L), function(a) {
        idx <- arm_subset(a)
        mean(dat[[col]][idx])
      }, numeric(1))

      obs_p0 <- vapply(0:(K - 1L), function(a) {
        idx <- arm_subset(a)
        mean(dat[[col]][idx] == 0)
      }, numeric(1))

      do.call(rbind, lapply(0:(K - 1L), function(a) {
        data.frame(
          endpoint = col,
          arm = a,
          input_baseline_mean = spec$baseline_mean,
          input_trt_logRR     = inp_eff[a + 1L],
          input_trt_mean      = inp_mean_vec[a + 1L],
          input_size          = spec$size,
          input_p_zero        = (spec$p_zero %||% 0),
          est_baseline_mean   = exp(b0),
          est_trt_logRR       = est_eff[a + 1L],
          est_size            = unname(fit$theta),
          obs_mean            = obs_mean_by_arm[a + 1L],
          obs_p0              = obs_p0[a + 1L],
          stringsAsFactors = FALSE
        )
      }))
    }))
    rownames(cnt_tbl) <- NULL
  }

  tte_tbl <- NULL
  if (length(idx_tte) > 0) {
    if (!requireNamespace("survival", quietly = TRUE)) {
      stop("survival is required for TTE endpoint summary via survival::coxph().")
    }

    tte_tbl <- do.call(rbind, lapply(seq_along(idx_tte), function(k) {
      j    <- idx_tte[k]
      col  <- endpoint_names[j]
      spec <- endpoint_details[[j]]

      censor_col <- paste0("Status_", k)
      if (!censor_col %in% names(dat)) stop("Missing censoring column ", censor_col, " for ", col, ".")

      time   <- dat[[col]]
      status <- dat[[censor_col]]

      fit <- if (has_trt && K > 1L) {
        survival::coxph(survival::Surv(time, status) ~ trtF)
      } else {
        survival::coxph(survival::Surv(time, status) ~ 1)
      }

      cf <- coef(fit)

      est_eff <- rep(0, K)
      names(est_eff) <- 0:(K - 1L)
      if (has_trt && K > 1L) {
        for (a in 1:(K - 1L)) {
          nm <- paste0("trtF", a)
          if (nm %in% names(cf)) est_eff[a + 1L] <- unname(cf[nm])
        }
      }

      inp_eff <- expand_input_effect(spec$trt_effect %||% NULL, K)

      obs_event_rate <- vapply(0:(K - 1L), function(a) {
        idx <- arm_subset(a)
        mean(status[idx])
      }, numeric(1))

      exp_rate <- vapply(0:(K - 1L), function(a) {
        idx <- arm_subset(a)
        sa <- status[idx]; ta <- time[idx]
        sum(sa) / sum(ta)
      }, numeric(1))

      do.call(rbind, lapply(0:(K - 1L), function(a) {
        data.frame(
          endpoint = col,
          arm = a,
          censor_col = censor_col,
          input_baseline_rate = spec$baseline_rate,
          input_trt_logHR     = inp_eff[a + 1L],
          input_trt_HR        = exp(inp_eff[a + 1L]), # new
          est_trt_logHR       = est_eff[a + 1L],
          est_trt_HR          = exp(est_eff[a + 1L]), # new
          obs_event_rate      = obs_event_rate[a + 1L],
          exp_rate            = exp_rate[a + 1L],
          stringsAsFactors = FALSE
        )
      }))
    }))
    rownames(tte_tbl) <- NULL
  }

  out <- list(
    target_correlation = cor_target,
    estimated_correlation_by_arm = cor_by_arm,
    continuous = cont_tbl,
    binary     = bin_tbl,
    count      = cnt_tbl,
    tte        = tte_tbl,
    n_arms = K
  )
  class(out) <- "summary.makeDataSim"
  out
}


#' @exportS3Method
print.summary.makeDataSim <- function(x, ...) {
  cat("<summary.makeDataSim>\n\n")
  cat("n_arms:", x$n_arms, "\n\n")

  cat("Target correlation:\n")
  print(x$target_correlation)

  cat("\nEstimated correlation (by arm):\n")
  for (nm in names(x$estimated_correlation_by_arm)) {
    cat("\n", nm, ":\n", sep = "")
    print(x$estimated_correlation_by_arm[[nm]])
  }

  if (!is.null(x$continuous)) { cat("\nContinuous endpoints (endpoint x arm):\n"); print(x$continuous) }
  if (!is.null(x$binary))     { cat("\nBinary endpoints (endpoint x arm):\n");     print(x$binary) }
  if (!is.null(x$count))      { cat("\nCount endpoints (endpoint x arm):\n");      print(x$count) }
  if (!is.null(x$tte))        { cat("\nTTE endpoints (endpoint x arm):\n");        print(x$tte) }

  invisible(x)
}


#' Pairwise visualization of simulated endpoints from a \code{makeDataSim}
#' object
#' 
#' \code{plot.makeDataSim()} produces a pairwise plot matrix for the simulated
#' endpoints stored in a \code{"makeDataSim"} object returned by
#' \code{\link{makeData}}.
#' 
#' The method uses \code{GGally::ggpairs()} to display scatterplots in the
#' upper triangle and pairwise correlations in the lower triangle for the
#' selected study arm. This is intended as a quick diagnostic tool for
#' inspecting the joint structure of simulated endpoints.
#' 
#' 
#' @param x An object of class \code{"makeDataSim"}, created by
#' \code{\link{makeData}}.
#' @param arm Numeric scalar indicating which treatment arm to plot.
#' 
#' If the simulated dataset includes a \code{trt} column, \code{arm} must be
#' one of \code{0, 1, \dots, K - 1}, where \code{K} is the number of arms.
#' 
#' If the object was generated in control-only mode (i.e., no \code{trt} column
#' is present), the \code{arm} argument is ignored and all simulated
#' observations are plotted.
#' @param names Optional character vector used to relabel the plotted endpoint
#' columns.
#' 
#' If supplied, \code{names} must have length equal to the number of simulated
#' endpoints. This affects only the displayed labels in the plot and does not
#' modify the underlying data.
#' @param \dots Additional arguments passed directly to
#' \code{GGally::ggpairs()}. This can be used to customize panels, labels,
#' sizing, or other plot options.
#' @return A \code{GGally::ggpairs} object, which is also a
#' \code{ggplot}-compatible object and can be printed or further modified using
#' standard \pkg{ggplot2} syntax.
#' @section Displayed data:
#' 
#' The plotting method extracts the simulated endpoint columns recorded in
#' \code{x$meta$endpoint_names}. These are the automatically renamed endpoint
#' variables, such as: \describe{ \item{Continuous endpoints}{\code{Cont_1},
#' \code{Cont_2}, \dots} \item{Binary endpoints}{\code{Bin_1}, \code{Bin_2},
#' \dots} \item{Count endpoints}{\code{Int_1}, \code{Int_2}, \dots}
#' \item{Time-to-event endpoints}{\code{TTE_1}, \code{TTE_2}, \dots} } Columns
#' such as \code{trt}, \code{Status_*}, and \code{enrollTime} are not included
#' in the pair plot.
#' @section Panel layout: The returned \code{ggpairs} object is configured as
#' follows: \itemize{ \item upper triangle: scatterplots for continuous-style
#' panels, \item lower triangle: pairwise correlations, \item title: \code{"Arm
#' <k>"} for the selected arm. } A \code{ggplot2::theme_bw()} theme is applied
#' by default.
#' @section Dependencies: This method requires both \pkg{GGally} and
#' \pkg{ggplot2}. If either package is not installed, the method throws an
#' error with an installation message.
#' @seealso \code{\link[=summary.makeDataSim]{summary.makeDataSim}} for
#' numerical summaries of the simulated data.
#' @examples
#' 
#' 
#' ep1 <- list(
#'   endpoint_type = "continuous",
#'   baseline_mean = 10,
#'   sd            = 2,
#'   trt_effect    = -1
#' )
#' 
#' ep2 <- list(
#'   endpoint_type = "binary",
#'   baseline_prob = 0.30,
#'   trt_prob      = 0.45
#' )
#' 
#' ep3 <- list(
#'   endpoint_type = "count",
#'   baseline_mean = 8,
#'   trt_count     = 10,
#'   size          = 20
#' )
#' 
#' R3 <- corr_make(
#'   num_endpoints = 3,
#'   values = rbind(
#'     c(1, 2, 0.20),
#'     c(1, 3, 0.10),
#'     c(2, 3, 0.15)
#'   )
#' )
#' 
#' sim_obj <- makeData(
#'   correlation_matrix    = R3,
#'   sample_size_per_group = 500,
#'   SEED                  = 1,
#'   endpoint_details      = list(ep1, ep2, ep3)
#' )
#' 
#' ## Plot control arm
#' plot(sim_obj)
#' 
#' ## Plot treatment arm with custom labels
#' plot(sim_obj, arm = 1, names = c("Biomarker", "Responder", "Hospitalizations"))
#' 
#' @exportS3Method
plot.makeDataSim <- function(x,
                             arm = 0,
                             names = NULL,
                             ...) {

  if (!requireNamespace("GGally", quietly = TRUE)) {
    stop("Package 'GGally' is required for plot.makeDataSim(). Please install it.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot.makeDataSim(). Please install it.")
  }

  dat  <- x$data
  meta <- x$meta

  if (is.null(meta$endpoint_names)) {
    stop("Object metadata does not include `endpoint_names`. Cannot determine endpoint columns.")
  }

  has_trt <- "trt" %in% names(dat)
  K <- meta$n_arms %||% if (has_trt) length(unique(dat$trt)) else 1L

  if (!has_trt) {
    dat_sub <- dat[, meta$endpoint_names, drop = FALSE]
    arm_lbl <- 0L
  } else {
    if (!is.numeric(arm) || length(arm) != 1 || is.na(arm)) stop("`arm` must be a single number.")
    arm <- as.integer(arm)
    if (arm < 0 || arm > (K - 1L)) stop("`arm` must be in {0,1,...,", (K - 1L), "}.")
    dat_sub <- dat[dat$trt == arm, meta$endpoint_names, drop = FALSE]
    arm_lbl <- arm
  }

  if (!is.null(names)) {
    if (!is.character(names)) stop("`names` must be a character vector.")
    if (length(names) != ncol(dat_sub)) {
      stop("`names` must have length equal to the number of endpoints (", ncol(dat_sub), ").")
    }
    colnames(dat_sub) <- names
  }

  p <- GGally::ggpairs(
    dat_sub,
    upper = list(continuous = "points"),
    lower = list(continuous = "cor"),
    progress = FALSE,
    ...
  ) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0("Arm ", arm_lbl))

  p
}


# -----------------------------------------------------------------------------
# Argument checks for makeData()
# -----------------------------------------------------------------------------

check_makeData_args <- function(correlation_matrix,
                                SEED,
                                sample_size_per_group,
                                endpoint_details,
                                enrollment_details,
                                non_fatal_censors_fatal,
                                target_correlation,
                                calibration_control,
                                arm_mode = c("auto", "full", "control")) {

  arm_mode <- match.arg(arm_mode)

  # ---- endpoint_details ----------------------------------------------------
  if (!is.list(endpoint_details) || length(endpoint_details) < 1) {
    stop("Error: `endpoint_details` must be a non-empty list of endpoint spec lists.")
  }
  p <- length(endpoint_details)

  endpoint_types <- vapply(endpoint_details, function(e) {
    if (is.null(e$endpoint_type)) stop("Error: Each endpoint spec must include `endpoint_type`.")
    normalize_endpoint_type(e$endpoint_type)
  }, character(1))

  # ---- SEED ----------------------------------------------------------------
  if (!is.null(SEED)) {
    if (!is.numeric(SEED) || length(SEED) != 1 || is.na(SEED)) {
      stop("Error: `SEED` must be NULL or a single numeric value.")
    }
  }

  # ---- single-endpoint mode: correlation_matrix=NULL -----------------------
  single_endpoint_mode <- is.null(correlation_matrix)
  if (single_endpoint_mode) {
    if (p != 1L) {
      stop("Error: If `correlation_matrix=NULL`, you must supply exactly one endpoint (length(endpoint_details)=1).")
    }
  } else {
    if (!is.matrix(correlation_matrix) || nrow(correlation_matrix) != p || ncol(correlation_matrix) != p) {
      stop("Error: Dimensions of `correlation_matrix` must match the number of endpoints (", p, "x", p, ").")
    }
    if (max(abs(correlation_matrix - t(correlation_matrix))) > 1e-10) {
      stop("Error: `correlation_matrix` must be symmetric.")
    }
    if (max(abs(diag(correlation_matrix) - 1)) > 1e-10) {
      stop("Error: `correlation_matrix` must have 1s on the diagonal.")
    }
    if (any(correlation_matrix < -1 | correlation_matrix > 1, na.rm = TRUE)) {
      stop("Error: all entries of `correlation_matrix` must be in [-1, 1].")
    }
  }

  # ---- infer K from arm-specific args --------------------------------------
  active_len <- function(x, nm) {
    if (is.null(x)) return(0L)
    if (!is.numeric(x) || anyNA(x)) stop("Error: `", nm, "` must be numeric with no NA.")
    as.integer(length(x))
  }

  lens_effect <- vapply(endpoint_details, function(e) active_len(e$trt_effect %||% NULL, "trt_effect"), integer(1))
  lens_bprob  <- vapply(endpoint_details, function(e) active_len(e$trt_prob  %||% NULL, "trt_prob"),  integer(1))
  lens_cmean  <- vapply(endpoint_details, function(e) active_len(e$trt_count %||% NULL, "trt_count"), integer(1))

  any_arm_info <- any(c(lens_effect, lens_bprob, lens_cmean) > 0L)

  control_only <- switch(
    arm_mode,
    control = TRUE,
    full    = FALSE,
    auto    = !any_arm_info
  )

  if (control_only) {
    K <- 1L
  } else {
    max_len <- max(c(lens_effect, lens_bprob, lens_cmean))
    K <- max(2L, max_len + 1L)
  }

  if (!control_only) {
    validate_len <- function(len, nm, j) {
      if (len == 0L) return(invisible(TRUE))
      if (!(len %in% c(1L, K - 1L))) {
        stop("Error: endpoint ", j, " `", nm, "` must have length 1 or K-1 (K=", K, ").")
      }
      invisible(TRUE)
    }
    for (j in seq_along(endpoint_details)) {
      validate_len(lens_effect[j], "trt_effect", j)
      validate_len(lens_bprob[j],  "trt_prob",  j)
      validate_len(lens_cmean[j],  "trt_count", j)
    }
  } else {
    if (any(c(lens_effect, lens_bprob, lens_cmean) > 0L)) {
      stop("Error: control-only mode prohibits specifying `trt_effect`, `trt_prob`, or `trt_count`.")
    }
  }

  # ---- sample_size_per_group ----------------------------------------------
  if (!is.numeric(sample_size_per_group) || anyNA(sample_size_per_group)) {
    stop("Error: `sample_size_per_group` must be numeric with no NA.")
  }
  if (length(sample_size_per_group) == 1L) {
    if (sample_size_per_group <= 0) stop("Error: `sample_size_per_group` must be > 0.")
    n_by_arm <- rep(as.integer(sample_size_per_group), K)
  } else if (length(sample_size_per_group) == K) {
    if (any(sample_size_per_group <= 0)) stop("Error: all `sample_size_per_group` entries must be > 0.")
    n_by_arm <- as.integer(sample_size_per_group)
  } else {
    stop("Error: `sample_size_per_group` must have length 1 or length K (K=", K, ").")
  }

  # ---- per-endpoint checks -------------------------------------------------
  for (j in seq_along(endpoint_details)) {
    spec <- endpoint_details[[j]]
    typ  <- endpoint_types[j]

    if (typ == "binary") {
      if (!is.null(spec$trt_effect) && !is.null(spec$trt_prob)) {
        stop("Error: Binary endpoint ", j, ": specify only one of `trt_effect` or `trt_prob`.")
      }
    }
    if (typ == "count") {
      if (!is.null(spec$trt_effect) && !is.null(spec$trt_count)) {
        stop("Error: Count endpoint ", j, ": specify only one of `trt_effect` or `trt_count`.")
      }
    }

    if (typ == "continuous") {
      if (is.null(spec$baseline_mean) || is.null(spec$sd)) {
        stop("Error: Continuous endpoint ", j, " requires `baseline_mean` and `sd`.")
      }
      if (!is.numeric(spec$baseline_mean) || length(spec$baseline_mean) != 1 || is.na(spec$baseline_mean)) {
        stop("Error: Continuous endpoint ", j, " `baseline_mean` must be a single numeric value.")
      }
      if (!is.numeric(spec$sd) || anyNA(spec$sd)) {
        stop("Error: Continuous endpoint ", j, " `sd` must be numeric with no NA.")
      }
      if (!(length(spec$sd) %in% c(1L, K))) {
        stop("Error: Continuous endpoint ", j, " `sd` must have length 1 or K (K=", K, ").")
      }
      if (any(spec$sd <= 0)) {
        stop("Error: Continuous endpoint ", j, " `sd` must be > 0.")
      }

    } else if (typ == "binary") {
      if (is.null(spec$baseline_prob)) stop("Error: Binary endpoint ", j, " requires `baseline_prob`.")
      if (!is.numeric(spec$baseline_prob) || length(spec$baseline_prob) != 1 || is.na(spec$baseline_prob)) {
        stop("Error: Binary endpoint ", j, " `baseline_prob` must be a single numeric value.")
      }
      if (!(spec$baseline_prob > 0 && spec$baseline_prob < 1)) {
        stop("Error: Binary endpoint ", j, " `baseline_prob` must be in (0, 1).")
      }
      if (!is.null(spec$trt_prob)) {
        if (!is.numeric(spec$trt_prob) || anyNA(spec$trt_prob)) stop("Error: Binary endpoint ", j, " `trt_prob` must be numeric with no NA.")
        if (!(length(spec$trt_prob) %in% c(1L, K - 1L))) stop("Error: Binary endpoint ", j, " `trt_prob` must have length 1 or K-1 (K=", K, ").")
        if (any(spec$trt_prob <= 0 | spec$trt_prob >= 1)) stop("Error: Binary endpoint ", j, " `trt_prob` must be in (0,1).")
      }

    } else if (typ == "count") {
      if (is.null(spec$baseline_mean) || is.null(spec$size)) {
        stop("Error: Count endpoint ", j, " requires `baseline_mean` and `size`.")
      }
      if (!is.numeric(spec$baseline_mean) || length(spec$baseline_mean) != 1 ||
          is.na(spec$baseline_mean) || spec$baseline_mean <= 0) {
        stop("Error: Count endpoint ", j, " `baseline_mean` must be a single positive numeric value.")
      }
      if (!is.numeric(spec$size) || length(spec$size) != 1 || is.na(spec$size) || spec$size <= 0) {
        stop("Error: Count endpoint ", j, " `size` must be a single positive numeric value.")
      }
      if (!is.null(spec$p_zero)) {
        if (!is.numeric(spec$p_zero) || length(spec$p_zero) != 1 || is.na(spec$p_zero)) {
          stop("Error: Count endpoint ", j, " `p_zero` must be a single numeric value.")
        }
        if (spec$p_zero < 0 || spec$p_zero > 1) {
          stop("Error: Count endpoint ", j, " `p_zero` must be in [0, 1].")
        }
      }
      if (!is.null(spec$trt_count)) {
        if (!is.numeric(spec$trt_count) || anyNA(spec$trt_count)) stop("Error: Count endpoint ", j, " `trt_count` must be numeric with no NA.")
        if (!(length(spec$trt_count) %in% c(1L, K - 1L))) stop("Error: Count endpoint ", j, " `trt_count` must have length 1 or K-1 (K=", K, ").")
        if (any(spec$trt_count <= 0)) stop("Error: Count endpoint ", j, " `trt_count` must be > 0.")
      }

    } else if (typ == "time-to-event") {
      if (is.null(spec$baseline_rate)) stop("Error: TTE endpoint ", j, " requires `baseline_rate`.")
      if (!is.numeric(spec$baseline_rate) || length(spec$baseline_rate) != 1 ||
          is.na(spec$baseline_rate) || spec$baseline_rate <= 0) {
        stop("Error: TTE endpoint ", j, " `baseline_rate` must be a single positive numeric value.")
      }
      if (!is.null(spec$censoring_rate)) {
        if (!is.numeric(spec$censoring_rate) || length(spec$censoring_rate) != 1 ||
            is.na(spec$censoring_rate) || spec$censoring_rate < 0) {
          stop("Error: TTE endpoint ", j, " `censoring_rate` must be a single numeric value >= 0.")
        }
      }
      if (!is.null(spec$fatal_event)) {
        if (!is.logical(spec$fatal_event) || length(spec$fatal_event) != 1) {
          stop("Error: TTE endpoint ", j, " `fatal_event` must be TRUE/FALSE.")
        }
      }
    }
  }

  # ---- fatal event structure checks ----------------------------------------
  tte_idx <- which(endpoint_types == "time-to-event")
  if (length(tte_idx) > 0) {
    fatal_flags <- vapply(tte_idx, function(j) isTRUE(endpoint_details[[j]]$fatal_event %||% FALSE), logical(1))

    if (sum(fatal_flags) > 2) stop("Error: Only two terminal (fatal) events are allowed.")

    if (sum(fatal_flags) > 0) {
      fatal_pos <- which(fatal_flags)
      if (any(fatal_pos > 2)) stop("Error: Terminal events must be the first or first and second TTE endpoints.")
    }

    if (any(diff(as.integer(fatal_flags)) == 1L)) {
      stop("Error: Fatal TTE endpoints must be listed before any non-fatal TTE endpoints (among TTE endpoints).")
    }
  }

  # ---- non_fatal_censors_fatal ---------------------------------------------
  if (!is.logical(non_fatal_censors_fatal) || length(non_fatal_censors_fatal) != 1) {
    stop("Error: `non_fatal_censors_fatal` must be TRUE or FALSE.")
  }

  # ---- enrollment / admin censoring ----------------------------------------
  if (!is.list(enrollment_details)) stop("Error: `enrollment_details` must be a list.")

  adm <- enrollment_details$administrative_censoring %||% NULL
  dist <- enrollment_details$enrollment_distribution %||% "none"

  valid_dists <- c("none","uniform","exponential","piecewise")
  if (!dist %in% valid_dists) {
    stop("Error: `enrollment_distribution` must be one of: ", paste(valid_dists, collapse = ", "))
  }

  if (!is.null(adm)) {
    if (!is.numeric(adm) || length(adm) != 1 || is.na(adm)) stop("Error: `administrative_censoring` must be numeric.")
    if (adm < 0) stop("Error: `administrative_censoring` must be >= 0.")
  }

  if (dist != "none" && is.null(adm)) {
    stop("Error: Please set `administrative_censoring` when enrollment_distribution != 'none'.")
  }

  if (dist == "exponential") {
    rate <- enrollment_details$enrollment_exponential_rate %||% NULL
    if (is.null(rate)) stop("Error: `enrollment_exponential_rate` must be specified when enrollment_distribution='exponential'.")
    if (!is.numeric(rate) || length(rate) != 1 || is.na(rate) || rate <= 0) {
      stop("Error: `enrollment_exponential_rate` must be a single positive number.")
    }
  }

  if (dist == "piecewise") {
    cuts  <- enrollment_details$piecewise_enrollment_cutpoints %||% NULL
    rates <- enrollment_details$piecewise_enrollment_rates %||% NULL

    if (is.null(cuts) || is.null(rates)) {
      stop("Error: `piecewise_enrollment_cutpoints` and `piecewise_enrollment_rates` must both be specified when enrollment_distribution='piecewise'.")
    }
    if (!is.numeric(cuts) || !is.numeric(rates)) {
      stop("Error: `piecewise_enrollment_cutpoints` and `piecewise_enrollment_rates` must be numeric vectors.")
    }
    if (any(diff(cuts) <= 0)) stop("Error: `piecewise_enrollment_cutpoints` must be strictly increasing.")
    if (length(rates) != length(cuts) - 1) stop("Error: length(piecewise_enrollment_rates) must be length(piecewise_enrollment_cutpoints) - 1.")
    if (any(rates <= 0)) stop("Error: all `piecewise_enrollment_rates` must be > 0.")
  }

  # ---- calibration checks --------------------------------------------------
  if (!is.logical(target_correlation) || length(target_correlation) != 1) {
    stop("Error: `target_correlation` must be TRUE/FALSE.")
  }

  if (isTRUE(target_correlation) && !single_endpoint_mode) {
    if (!is.list(calibration_control)) stop("Error: `calibration_control` must be a list when target_correlation=TRUE.")

    # New names (preferred) + backward-compatible old names
    n_mc      <- calibration_control$n_mc      %||% calibration_control$n_obs      %||% 10000
    tol       <- calibration_control$tol       %||% 0.001
    maxit     <- calibration_control$maxit     %||% 100
    rho_cap   <- calibration_control$rho_cap   %||% calibration_control$rho_max    %||% 0.999
    ensure_pd <- calibration_control$ensure_pd %||% calibration_control$ensure_cor_mat %||% TRUE
    convT     <- calibration_control$conv_norm_type %||% "F"

    if (!is.numeric(n_mc) || length(n_mc) != 1 || is.na(n_mc) || n_mc <= 0)
      stop("Error: `calibration_control$n_mc` (or legacy `n_obs`) must be > 0.")
    if (!is.numeric(tol) || length(tol) != 1 || is.na(tol) || tol <= 0)
      stop("Error: `calibration_control$tol` must be > 0.")
    if (!is.numeric(maxit) || length(maxit) != 1 || is.na(maxit) || maxit <= 0)
      stop("Error: `calibration_control$maxit` must be > 0.")
    if (!is.numeric(rho_cap) || length(rho_cap) != 1 || is.na(rho_cap) || rho_cap <= 0 || rho_cap >= 1) {
      stop("Error: `calibration_control$rho_cap` (or legacy `rho_max`) must be in (0, 1).")
    }
    if (!is.logical(ensure_pd) || length(ensure_pd) != 1)
      stop("Error: `calibration_control$ensure_pd` (or legacy `ensure_cor_mat`) must be TRUE/FALSE.")
    if (!is.character(convT) || length(convT) != 1)
      stop("Error: `calibration_control$conv_norm_type` must be a single character value.")
  }

  invisible(list(
    n_arms = K,
    n_by_arm = n_by_arm,
    endpoint_types = endpoint_types,
    control_only = control_only,
    single_endpoint_mode = single_endpoint_mode
  ))
}
