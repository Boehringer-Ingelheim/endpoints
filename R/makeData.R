#' Simulate trial data with mixed endpoint types, optional Gaussian-copula
#' dependence, and enrollment/censoring features
#'
#' \code{makeData()} simulates subject-level trial data with one or more
#' endpoints. Supported endpoint types include continuous, binary, count, and
#' time-to-event outcomes. When multiple endpoints are supplied, dependence may
#' be induced through a Gaussian copula (NORTA-style construction), with
#' optional numerical calibration to approximately match a target Pearson
#' correlation matrix on the observed endpoint scale.
#'
#' The function also supports:
#' \itemize{
#'   \item multiple treatment arms,
#'   \item independent censoring for time-to-event outcomes,
#'   \item fatal and non-fatal time-to-event logic (including semi-competing
#'   risks),
#'   \item administrative censoring,
#'   \item stochastic enrollment (uniform, exponential, or piecewise
#'   exponential), and
#'   \item generation of longitudinal data.
#' }
#'
#' @param correlation_matrix A numeric correlation matrix specifying the
#'   dependence structure across endpoints. If \code{target_correlation = FALSE},
#'   this is interpreted on the \emph{latent Gaussian} scale. If \code{NULL},
#'   \code{makeData()} enters single-endpoint mode, and exactly one endpoint
#'   must be supplied in \code{endpoint_details}. Typically created with
#'   \code{\link{corr_make}}.
#'
#' @param SEED Optional numeric scalar seed used to initialize the
#'   random-number generator via \code{set.seed()}. If \code{NULL}, the current
#'   RNG state is used.
#'
#' @param sample_size_per_group Integer scalar or integer vector giving the
#'   sample size per arm. If a scalar, the same sample size is used for all
#'   arms. If a vector, it must be of length \eqn{K}, where \eqn{K} is the
#'   total number of arms, including control.
#'
#' @param endpoint_details A non-empty list of endpoint specification lists,
#'   one per endpoint.
#'
#'   Each sub-list defines the marginal distribution and arm-specific
#'   parameters for one endpoint. See \code{\link{endpoint_details}} for the
#'   complete schema, supported endpoint types, and detailed examples.
#'
#' @param enrollment_details A named list controlling enrollment and
#'   administrative follow-up. Fields include:
#'   \describe{
#'     \item{administrative_censoring}{
#'       Numeric scalar. If non-\code{NULL}, follow-up is administratively
#'       censored at this time.
#'     }
#'     \item{enrollment_distribution}{
#'       Character. One of \code{"none"}, \code{"uniform"},
#'       \code{"exponential"}, or \code{"piecewise"}.
#'     }
#'     \item{enrollment_exponential_rate}{
#'       Numeric scalar rate for exponential enrollment when
#'       \code{enrollment_distribution = "exponential"}.
#'     }
#'     \item{piecewise_enrollment_cutpoints}{
#'       Numeric vector of cutpoints defining intervals for piecewise
#'       enrollment.
#'     }
#'     \item{piecewise_enrollment_rates}{
#'       Numeric vector of rates, one per interval, for piecewise exponential
#'       enrollment.
#'     }
#'   }
#'
#'   See \code{\link{enrollment_details}} for full details and recommended
#'   parameterization.
#'
#' @param non_fatal_censors_fatal Logical. Controls semi-competing risks
#'   behavior when multiple TTE endpoints are present. If \code{TRUE}, censoring
#'   of a non-fatal endpoint can censor subsequent TTE endpoints according to
#'   the package's semi-competing-risk rules. If \code{FALSE}, non-fatal
#'   censoring does not censor other TTE endpoints.
#'
#' @param target_correlation Logical scalar.
#'
#'   If \code{TRUE} and \code{correlation_matrix} is not \code{NULL}, a
#'   numerical calibration step is used to search for a latent Gaussian
#'   correlation matrix whose transformed endpoints approximately match the
#'   requested observed Pearson correlation matrix.
#'
#'   If \code{FALSE}, the supplied \code{correlation_matrix} is used directly
#'   as the latent Gaussian correlation matrix.
#'
#'   This value should generally be set to \code{TRUE}, unless there is a
#'   compelling reason otherwise.
#'
#' @param arm_mode Character string controlling how the number of treatment
#'   arms is determined. Must be one of:
#'   \describe{
#'     \item{auto}{
#'       Infer the number of arms from the endpoint specifications. If no
#'       treatment-specific quantities are supplied, control-only mode is used.
#'     }
#'     \item{full}{
#'       Force a full multi-arm interpretation, even if some endpoint
#'       specifications omit treatment-specific effects.
#'     }
#'     \item{control}{
#'       Force control-only mode. In this case, treatment-specific arguments
#'       such as \code{trt_effect}, \code{trt_prob}, and \code{trt_count} are
#'       not allowed.
#'     }
#'   }
#'
#'   This value should generally be set to \code{"auto"}.
#'
#' @param calibration_control Named list of controls for the correlation
#'   calibration routine used when \code{target_correlation = TRUE}. Common
#'   controls include Monte Carlo size, tolerance, iteration caps, and
#'   constraints to ensure a valid correlation matrix. See
#'   \code{\link{calibration_control}} for more detail.
#'
#' @return An object of class \code{"makeDataSim"}.
#'
#'   Internally, this is a list with at least:
#'   \describe{
#'     \item{data}{
#'       A \code{data.frame} containing the simulated dataset.
#'     }
#'     \item{meta}{
#'       A metadata list storing endpoint types, endpoint names, arm counts,
#'       and input settings.
#'     }
#'   }
#'
#'   The \code{data} component may include:
#'   \describe{
#'     \item{trt}{
#'       Treatment arm indicator, omitted in control-only mode.
#'     }
#'     \item{Cont_1, Cont_2, \dots}{
#'       Continuous endpoints.
#'     }
#'     \item{Bin_1, Bin_2, \dots}{
#'       Binary endpoints.
#'     }
#'     \item{Int_1, Int_2, \dots}{
#'       Count endpoints.
#'     }
#'     \item{TTE_1, TTE_2, \dots}{
#'       Observed time-to-event variables.
#'     }
#'     \item{Status_1, Status_2, \dots}{
#'       TTE event indicators, with \code{1 = event} and
#'       \code{0 = censored}.
#'     }
#'     \item{enrollTime}{
#'       Enrollment times, when stochastic enrollment is used.
#'     }
#'   }
#'
#' @section Overview:
#' When multiple endpoints are simulated, \code{makeData()} uses a Gaussian
#' copula construction. If \eqn{\mathbf{Z} \sim N(\mathbf{0}, \Sigma_Z)} is a
#' latent multivariate normal vector, then each endpoint is generated as
#' \deqn{
#' X_j = F_j^{-1}\{\Phi(Z_j)\},
#' }
#' where \eqn{\Phi} is the standard normal CDF and \eqn{F_j^{-1}} is the
#' endpoint-specific quantile function.
#'
#' For continuous endpoints this yields Gaussian margins; for binary, count,
#' and time-to-event endpoints, the corresponding marginal quantile functions
#' are applied to the copula uniforms.
#'
#' @section Treatment arms:
#' The total number of study arms is generally determined from the lengths of
#' treatment-specific inputs in \code{endpoint_details}, for example the length
#' of \code{trt_effect}. When treatment arms are present, the output includes a
#' \code{trt} column coded as \code{0, 1, 2, \dots}. This can also be
#' controlled via \code{arm_mode}.
#'
#' @section Administrative censoring and enrollment:
#' If \code{administrative_censoring} is supplied, all TTE outcomes are
#' truncated at the maximum available follow-up. If stochastic enrollment is
#' enabled, each subject receives an \code{enrollTime}, and maximum observable
#' follow-up is reduced to \eqn{\mathcal{A} - T_E}, where \eqn{\mathcal{A}} is
#' the administrative censoring time and \eqn{T_E} is the enrollment time.
#'
#' @section Output object:
#' The returned object has class \code{"makeDataSim"} and is intended to be
#' used with:
#' \itemize{
#'   \item \code{print()} for a compact overview,
#'   \item \code{summary()} for marginal and correlation diagnostics, and
#'   \item \code{plot()} for quick visualization of endpoint relationships.
#' }
#'
#' @seealso
#' \code{\link{endpoint_details}} for endpoint specification details.
#'
#' \code{\link{enrollment_details}} for administrative censoring and stochastic
#' enrollment options.
#'
#' \code{\link{calibration_control}} for calibration tuning parameters.
#'
#' \code{\link{corr_make}} for creating correlation matrices.
#'
#' See the \code{vignette("user_guide", package = "endpoints")} vignette for
#' introductory examples, and the
#' \code{vignette("longitudinal_data", package = "endpoints")} vignette for
#' simulating longitudinal data.
#'
#' @references
#' Cario, M. C., & Nelson, B. L. (1997). *Modeling and generating random
#' vectors with arbitrary marginal distributions and correlation matrix* (pp.
#' 1-19). Technical Report, Department of Industrial Engineering and Management
#' Sciences, Northwestern University, Evanston, Illinois.
#'
#' @keywords simulation
#'
#' @examples
#' library(endpoints)
#'
#' ## One continuous endpoint, two-arm trial
#' ep1 <- list(
#'   endpoint_type = "continuous",
#'   baseline_mean = 10,
#'   sd            = 2,
#'   trt_effect    = -1
#' )
#'
#' sim1 <- makeData(
#'   correlation_matrix    = NULL,
#'   sample_size_per_group = 200,
#'   SEED                  = 1,
#'   endpoint_details      = list(ep1)
#' )
#'
#' sim1
#' summary(sim1)
#'
#' ## Three correlated endpoints
#' R3 <- corr_make(
#'   num_endpoints = 3,
#'   values = rbind(
#'     c(1, 2, 0.20),
#'     c(1, 3, 0.10),
#'     c(2, 3, 0.15)
#'   )
#' )
#'
#' ep_cont <- list(
#'   endpoint_type = "continuous",
#'   baseline_mean = 10,
#'   sd            = 2,
#'   trt_effect    = -1
#' )
#'
#' ep_bin <- list(
#'   endpoint_type = "binary",
#'   baseline_prob = 0.30,
#'   trt_prob      = 0.45
#' )
#'
#' ep_cnt <- list(
#'   endpoint_type = "count",
#'   baseline_mean = 8,
#'   trt_count     = 10,
#'   size          = 20,
#'   p_zero        = 0
#' )
#'
#' sim3 <- makeData(
#'   correlation_matrix    = R3,
#'   sample_size_per_group = 1000,
#'   SEED                  = 123,
#'   endpoint_details      = list(ep_cont, ep_bin, ep_cnt),
#'   target_correlation    = TRUE
#' )
#'
#' summary(sim3)
#'
#' ## One TTE endpoint with administrative censoring and exponential enrollment
#' ep_tte <- list(
#'   endpoint_type  = "tte",
#'   baseline_rate  = 1 / 24,
#'   trt_effect     = log(0.8),
#'   fatal_event    = TRUE
#' )
#'
#' sim_tte <- makeData(
#'   correlation_matrix    = NULL,
#'   sample_size_per_group = 500,
#'   endpoint_details      = list(ep_tte),
#'   enrollment_details    = list(
#'     administrative_censoring    = 24,
#'     enrollment_distribution     = "exponential",
#'     enrollment_exponential_rate = 1 / 4
#'   )
#' )
#'
#' summary(sim_tte)
#'
#' @export
makeData <- function(
    correlation_matrix,
    SEED = NULL,
    sample_size_per_group,
    endpoint_details,
    enrollment_details = list(),
    non_fatal_censors_fatal = FALSE,
    target_correlation = TRUE,
    arm_mode = c("auto", "full", "control"),
    calibration_control = list(
      n_mc = 10000,
      tol = 0.001,
      maxit = 100,
      rho_cap = 0.999,
      ensure_pd = TRUE,
      conv_norm_type = "F"
    )
) {

  arm_mode <- match.arg(arm_mode)

  # ---- enrollment defaults -------------------------------------------------
  enrollment_details <- utils::modifyList(
    list(
      administrative_censoring       = NULL,
      enrollment_distribution        = "none",
      enrollment_exponential_rate    = NULL,
      piecewise_enrollment_cutpoints = NULL,
      piecewise_enrollment_rates     = NULL
    ),
    enrollment_details
  )
  enrollment_details$enrollment_distribution <- match.arg(
    enrollment_details$enrollment_distribution,
    c("none","uniform","exponential","piecewise")
  )

  # ---- checks (delegated to helper function) --------------------------------------------------
  chk <- check_makeData_args(
    correlation_matrix = correlation_matrix,
    SEED = SEED,
    sample_size_per_group = sample_size_per_group,
    endpoint_details = endpoint_details,
    enrollment_details = enrollment_details,
    non_fatal_censors_fatal = non_fatal_censors_fatal,
    target_correlation = target_correlation,
    calibration_control = calibration_control,
    arm_mode = arm_mode
  )

  endpoint_types <- chk$endpoint_types
  K              <- chk$n_arms
  n_by_arm       <- chk$n_by_arm
  control_only   <- isTRUE(chk$control_only)

  p <- length(endpoint_details)

  # ---- design skeleton: one row per arm, for the copula generation -----------------------------------
  tt <- data.frame(
    trt = 0:(K - 1L),
    n_i = as.integer(n_by_arm)
  )

  # ---- helpers -------------------------------------------------------------
  expand_eff_K <- function(eff, K) {
    if (is.null(eff)) return(rep(0, K))
    if (!is.numeric(eff) || anyNA(eff)) stop("`trt_effect` must be numeric (or NULL).")
    if (length(eff) == 1L) return(c(0, rep(eff, K - 1L)))
    if (length(eff) == (K - 1L)) return(c(0, eff))
    stop("`trt_effect` must have length 1 or K-1 (or be NULL).")
  }

  expand_param_K_including_control <- function(x, K, name = "param") {
    if (is.null(x)) stop("Internal error: missing ", name)
    if (!is.numeric(x) || anyNA(x)) stop("`", name, "` must be numeric with no NA.")
    if (length(x) == 1L) return(rep(x, K))
    if (length(x) == K)  return(as.numeric(x))
    stop("`", name, "` must have length 1 or K (K=", K, ").")
  }

  expand_active <- function(x, K, name = "arg") {
    if (is.null(x)) return(NULL)
    if (!is.numeric(x) || anyNA(x)) stop("`", name, "` must be numeric with no NA.")
    if (length(x) == 1L) return(rep(x, K - 1L))
    if (length(x) == (K - 1L)) return(as.numeric(x))
    stop("`", name, "` must have length 1 or K-1 (K=", K, ").")
  }

  # ---- build mu1..mu_p and (optionally) sd1..sd_p (for normal endpoints) --------------------------
  for (j in seq_len(p)) {
    spec <- endpoint_details[[j]]
    typ  <- endpoint_types[j]

    if (typ == "continuous") {
      eff_vec <- expand_eff_K(spec$trt_effect %||% NULL, K)
      arm_eff <- eff_vec[tt$trt + 1L]
      tt[[paste0("mu", j)]] <- spec$baseline_mean + arm_eff

      sd_vec <- expand_param_K_including_control(spec$sd, K, name = "sd")
      if (any(sd_vec <= 0)) stop("Continuous endpoint j=", j, " `sd` must be > 0.")
      tt[[paste0("sd", j)]] <- sd_vec[tt$trt + 1L]

    } else if (typ == "binary") {

      logit0 <- logit(spec$baseline_prob)

      if (!is.null(spec$trt_prob)) {
        p_trt <- expand_active(spec$trt_prob, K, name = "trt_prob")
        eff_active <- logit(p_trt) - logit0
        eff_vec <- c(0, eff_active)
      } else {
        eff_vec <- expand_eff_K(spec$trt_effect %||% NULL, K)
      }

      arm_eff <- eff_vec[tt$trt + 1L]
      tt[[paste0("mu", j)]] <- inv_logit(logit0 + arm_eff)

    } else if (typ == "count") {

      if (!is.null(spec$trt_count)) {
        mu_trt <- expand_active(spec$trt_count, K, name = "trt_count")
        eff_active <- log(mu_trt / spec$baseline_mean)
        eff_vec <- c(0, eff_active)
      } else {
        eff_vec <- expand_eff_K(spec$trt_effect %||% NULL, K)
      }

      arm_eff <- eff_vec[tt$trt + 1L]
      tt[[paste0("mu", j)]] <- spec$baseline_mean * exp(arm_eff)

    } else if (typ == "time-to-event") {

      eff_vec <- expand_eff_K(spec$trt_effect %||% NULL, K)
      arm_eff <- eff_vec[tt$trt + 1L]
      tt[[paste0("mu", j)]] <- spec$baseline_rate * exp(arm_eff)

    } else stop("Unknown endpoint type at j=", j)
  }

  if (!is.null(SEED)) set.seed(SEED)

  # ---- correlation handler -------------------------------------------------
  single_endpoint_mode <- is.null(correlation_matrix)

  if (!single_endpoint_mode) {
    L_latent_default <- make_pd(correlation_matrix)$L
  }

  # normalize calibration_control
  cc <- calibration_control
  cc_n_mc      <- cc$n_mc            %||% 10000
  cc_tol       <- cc$tol             %||% 0.001
  cc_maxit     <- cc$maxit           %||% 100
  cc_rho_cap   <- cc$rho_cap         %||% 0.999
  cc_ensure_pd <- cc$ensure_pd       %||% TRUE
  cc_conv_type <- cc$conv_norm_type  %||% "O"

  # ---- simulate per arm ----------------------------------------------------
  conditional_dist_sim <- function(i) {
    n_i <- tt$n_i[i]

    create_dist_function <- function(j, mu_value, sd_value = NULL) {
      spec <- endpoint_details[[j]]
      typ  <- endpoint_types[j]

      if (typ == "continuous") {
        local({
          m <- mu_value; s <- sd_value
          function(u) qnorm(u, mean = m, sd = s)
        })
      } else if (typ == "binary") {
        local({
          p <- mu_value
          function(u) qbinom(u, size = 1, prob = p)
        })
      } else if (typ == "count") {
        local({
          mu <- mu_value; size <- spec$size; p0 <- (spec$p_zero %||% 0)
          function(u) qzinb_mixture(u, mu = mu, size = size, p0 = p0)
        })
      } else if (typ == "time-to-event") {
        local({
          r <- mu_value
          function(u) qexp(u, rate = r)
        })
      } else stop("Unknown endpoint type at j=", j)
    }

    qfun_list <- lapply(seq_len(p), function(j) {
      mu_value <- tt[i, paste0("mu", j)]
      if (endpoint_types[j] == "continuous") {
        sd_value <- tt[i, paste0("sd", j)]
        create_dist_function(j, mu_value, sd_value = sd_value)
      } else {
        create_dist_function(j, mu_value)
      }
    })

    if (single_endpoint_mode) {
      U <- matrix(runif(n_i), ncol = 1)
      eps <- 1e-12
      U <- pmin(pmax(U, eps), 1 - eps)

    } else {
      if (isTRUE(target_correlation)) {
        latent_cor_cal <- calibrate_latent_cor_matrix(
          target_cor      = correlation_matrix,
          qfuns           = qfun_list,
          ensure_pd       = isTRUE(cc_ensure_pd),
          conv_norm_type  = cc_conv_type,
          return_diagnostics = FALSE,
          n_mc            = cc_n_mc,
          seed            = NULL,
          tol             = cc_tol,
          maxit           = cc_maxit,
          rho_cap         = cc_rho_cap
        )
        L_use <- make_pd(latent_cor_cal)$L
      } else {
        L_use <- L_latent_default
      }

      U <- draw_gaussian_copula_u(n = n_i, L = L_use, seed = NULL)
    }

    X <- lapply(seq_len(p), function(j) qfun_list[[j]](U[, j]))
    sim_data <- as.data.frame(X)
    names(sim_data) <- paste0("V", seq_len(p))

    if (!control_only) sim_data$trt <- tt[i, "trt"]
    sim_data
  }

  total_sim_data <- do.call(rbind, lapply(seq_len(nrow(tt)), conditional_dist_sim))

  # ---- TTE censoring + indicators -----------------------------------------
  tte_idx <- which(endpoint_types == "time-to-event")
  if (length(tte_idx) > 0) {

    censoring_rates <- vapply(tte_idx, function(j) endpoint_details[[j]]$censoring_rate %||% 0, numeric(1))
    fatal_events <- vapply(tte_idx, function(j) isTRUE(endpoint_details[[j]]$fatal_event %||% FALSE), logical(1))

    for (k in seq_along(tte_idx)) {
      j <- tte_idx[k]
      time_col <- paste0("V", j)

      ev_times <- total_sim_data[[time_col]]
      cr <- censoring_rates[k]
      cens_times <- if (cr <= 0) rep(Inf, length(ev_times)) else rexp(length(ev_times), rate = cr)

      is_censored <- ev_times > cens_times
      total_sim_data[[time_col]] <- pmin(ev_times, cens_times)
      total_sim_data[[paste0("Status_", k)]] <- as.integer(!is_censored)
    }

    if (non_fatal_censors_fatal) {
      fatal_k <- which(fatal_events)
      nonfatal_k <- setdiff(seq_along(tte_idx), fatal_k)

      for (k in nonfatal_k) {
        time_col_k <- paste0("V", tte_idx[k])
        cens_col_k <- paste0("Status_", k)

        is_cens <- total_sim_data[[cens_col_k]] == 0
        if (!any(is_cens)) next

        t_cens <- total_sim_data[[time_col_k]][is_cens]

        for (kk in seq_along(tte_idx)) {
          if (kk == k) next
          time_col_kk <- paste0("V", tte_idx[kk])
          cens_col_kk <- paste0("Status_", kk)

          total_sim_data[[time_col_kk]][is_cens] <- pmin(total_sim_data[[time_col_kk]][is_cens], t_cens)
          total_sim_data[[cens_col_kk]][is_cens] <- 0
        }
      }
    }

    if (length(tte_idx) > 1 && any(fatal_events)) {
      if (length(tte_idx) >= 2 && fatal_events[1] && fatal_events[2]) {
        t1 <- total_sim_data[[paste0("V", tte_idx[1])]]
        t2 <- total_sim_data[[paste0("V", tte_idx[2])]]
        c1 <- total_sim_data[["Status_1"]]
        c2 <- total_sim_data[["Status_2"]]

        cond1 <- (c1 == 1) & (t1 < t2)
        total_sim_data[[paste0("V", tte_idx[2])]][cond1] <- t1[cond1]
        total_sim_data[["Status_2"]][cond1] <- 0

        cond2 <- (c2 == 1) & (t2 < t1)
        total_sim_data[[paste0("V", tte_idx[1])]][cond2] <- t2[cond2]
        total_sim_data[["Status_1"]][cond2] <- 0

        cond3 <- (c1 == 0) & (c2 == 0)
        m <- pmin(t1[cond3], t2[cond3])
        total_sim_data[[paste0("V", tte_idx[1])]][cond3] <- m
        total_sim_data[[paste0("V", tte_idx[2])]][cond3] <- m
      }

      for (k in seq_along(tte_idx)) {
        if (!fatal_events[k]) next
        t_k <- total_sim_data[[paste0("V", tte_idx[k])]]

        if (k < length(tte_idx)) {
          for (kk in (k + 1):length(tte_idx)) {
            time_col_kk <- paste0("V", tte_idx[kk])
            cens_col_kk <- paste0("Status_", kk)

            t_kk <- total_sim_data[[time_col_kk]]
            cond <- t_kk > t_k

            total_sim_data[[time_col_kk]][cond] <- t_k[cond]
            total_sim_data[[cens_col_kk]][cond] <- 0
          }
        }
      }
    }
  }

  # ---- Enrollment + admin censoring ---------------------------------------
  admin_cens <- enrollment_details$administrative_censoring
  if (!is.null(admin_cens) && admin_cens > 0) {
    nSubs <- nrow(total_sim_data)
    enroll_time <- rep(0, nSubs)

    if (enrollment_details$enrollment_distribution == "uniform") {
      enroll_time <- runif(nSubs, min = 0, max = admin_cens)

    } else if (enrollment_details$enrollment_distribution == "exponential") {
      rate <- enrollment_details$enrollment_exponential_rate
      if (is.null(rate) || rate <= 0)
        stop("For exponential enrollment, provide `enrollment_exponential_rate > 0`.")
      tmp <- rexp(nSubs, rate = rate)
      enroll_time <- pmin(tmp, admin_cens)

    } else if (enrollment_details$enrollment_distribution == "piecewise") {
      cuts  <- enrollment_details$piecewise_enrollment_cutpoints
      rates <- enrollment_details$piecewise_enrollment_rates

      piecewise_draw <- function() {
        piece_len <- diff(cuts)
        t <- 0
        for (k in seq_along(piece_len)) {
          w <- rexp(1, rate = rates[k])
          if (w < piece_len[k]) { t <- t + w; break }
          t <- t + piece_len[k]
        }
        min(t, max(cuts))
      }
      enroll_time <- replicate(nSubs, piecewise_draw())
    }

    total_sim_data$enrollTime <- enroll_time

    if (length(tte_idx) > 0) {
      max_follow_up <- pmax(0, admin_cens - enroll_time)
      for (k in seq_along(tte_idx)) {
        time_col <- paste0("V", tte_idx[k])
        cens_col <- paste0("Status_", k)

        t <- total_sim_data[[time_col]]
        is_admin_cens <- t > max_follow_up

        total_sim_data[[time_col]] <- pmin(t, max_follow_up)
        total_sim_data[[cens_col]][is_admin_cens] <- 0
      }
    }
  }

  # ---- rename endpoints ----------------------------------------------------
  cont_k <- 0; bin_k <- 0; tte_k <- 0; cnt_k <- 0
  new_names <- character(p)

  for (j in seq_len(p)) {
    if (endpoint_types[j] == "continuous")    { cont_k <- cont_k + 1; new_names[j] <- paste0("Cont_", cont_k) }
    if (endpoint_types[j] == "binary")        { bin_k  <- bin_k  + 1; new_names[j] <- paste0("Bin_",  bin_k) }
    if (endpoint_types[j] == "count")         { cnt_k  <- cnt_k  + 1; new_names[j] <- paste0("Int_",  cnt_k) }
    if (endpoint_types[j] == "time-to-event") { tte_k  <- tte_k  + 1; new_names[j] <- paste0("TTE_",  tte_k) }
  }
  names(total_sim_data)[match(paste0("V", seq_len(p)), names(total_sim_data))] <- new_names

  meta <- list(
    correlation_matrix   = correlation_matrix,
    target_correlation   = isTRUE(target_correlation) && !single_endpoint_mode,
    endpoint_details     = endpoint_details,
    endpoint_types       = endpoint_types,
    endpoint_names       = new_names,
    n_arms               = K,
    n_by_arm             = as.integer(n_by_arm),
    control_only         = control_only,
    single_endpoint_mode = single_endpoint_mode
  )

  new_makeDataSim(total_sim_data, meta)
}
