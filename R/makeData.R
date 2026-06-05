#' Simulate trial data with mixed endpoint types, optional Gaussian-copula
#' dependence, stochastic enrollment, and trial-calendar censoring
#'
#' \code{makeData()} simulates subject-level trial data with one or more
#' endpoints. Supported endpoint types include continuous, binary, count, and
#' time-to-event outcomes. When multiple endpoints are supplied, dependence may
#' be induced through a Gaussian copula, with optional numerical calibration to
#' approximately match a target Pearson correlation matrix on the observed
#' endpoint scale.
#'
#' The function supports multiple treatment arms, independent censoring for
#' time-to-event outcomes, fatal and non-fatal time-to-event logic, stochastic
#' enrollment, maximum follow-up caps, variable-duration trial endings, fixed
#' calendar trial endings, and event-driven trial endings.
#'
#' @param correlation_matrix A numeric correlation matrix specifying the
#'   dependence structure across endpoints. If \code{target_correlation = FALSE},
#'   this is interpreted on the latent Gaussian scale. If \code{NULL},
#'   \code{makeData()} enters single-endpoint mode.
#'
#' @param SEED Optional numeric scalar seed used to initialize the random-number
#'   generator via \code{set.seed()}.
#'
#' @param sample_size_per_group Integer scalar or integer vector giving the
#'   sample size per arm.
#'
#' @param endpoint_details A non-empty list of endpoint specification lists, one
#'   per endpoint.
#'
#' @param enrollment_details A named list controlling enrollment. Fields include:
#'   \describe{
#'     \item{enrollment_distribution}{
#'       Character. One of \code{"none"}, \code{"exponential"}, or
#'       \code{"piecewise"}. The value \code{"exponential"} denotes homogeneous
#'       Poisson-process accrual, with exponentially distributed inter-arrival
#'       times. The value \code{"piecewise"} denotes piecewise homogeneous
#'       Poisson-process accrual.
#'     }
#'     \item{enrollment_exponential_rate}{
#'       Numeric scalar accrual rate per unit time when
#'       \code{enrollment_distribution = "exponential"}.
#'     }
#'     \item{piecewise_enrollment_cutpoints}{
#'       Numeric vector of cutpoints defining calendar intervals for piecewise
#'       enrollment. Should start at 0.
#'     }
#'     \item{piecewise_enrollment_rates}{
#'       Numeric vector of accrual rates, one per interval, for piecewise
#'       Poisson-process enrollment.
#'     }
#'   }
#'
#' @param followup_details A named list controlling subject-level follow-up
#'   limits. Fields include:
#'   \describe{
#'     \item{min_followup}{
#'       Planned minimum follow-up used by trial-ending rules such as
#'       \code{"last_patient_min_followup"}.
#'     }
#'     \item{max_followup}{
#'       Maximum observable follow-up for any individual subject.
#'     }
#'   }
#'
#' @param trial_end_details A named list controlling the calendar time of final
#'   analysis or database cutoff. Fields include:
#'   \describe{
#'     \item{type}{
#'       Character. One of \code{"none"}, \code{"fixed_calendar"},
#'       \code{"last_patient_min_followup"}, or \code{"event_driven"}.
#'     }
#'     \item{trial_end_time}{
#'       Numeric scalar calendar time for \code{type = "fixed_calendar"}.
#'     }
#'     \item{event_endpoint}{
#'       TTE endpoint used for event-driven stopping. May be a TTE ordinal such
#'       as \code{1}, or a character value such as \code{"TTE_1"}.
#'     }
#'     \item{target_events}{
#'       Positive integer number of events required for
#'       \code{type = "event_driven"}.
#'     }
#'     \item{require_min_followup}{
#'       Logical. For event-driven trials, if \code{TRUE}, the trial cannot end
#'       before the last randomized subject has \code{min_followup}.
#'     }
#'     \item{max_trial_duration}{
#'       Optional positive scalar maximum calendar trial duration.
#'     }
#'     \item{target_not_reached}{
#'       Character. One of \code{"error"}, \code{"max_trial_duration"}, or
#'       \code{"last_patient_min_followup"}.
#'     }
#'   }
#'
#' @param non_fatal_censors_fatal Logical. Controls semi-competing risks
#'   behavior when multiple TTE endpoints are present.
#'
#' @param target_correlation Logical scalar. If \code{TRUE}, numerical
#'   calibration is used to approximately match the requested observed Pearson
#'   correlation matrix.
#'
#' @param arm_mode Character string controlling how the number of treatment arms
#'   is determined. Must be one of \code{"auto"}, \code{"full"}, or
#'   \code{"control"}.
#'
#' @param calibration_control Named list of controls for the correlation
#'   calibration routine.
#'
#' @return An object of class \code{"makeDataSim"}.
#'
#' @export
makeData <- function(
    correlation_matrix,
    SEED = NULL,
    sample_size_per_group,
    endpoint_details,
    enrollment_details = list(),
    followup_details = list(),
    trial_end_details = list(),
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

  # ---- calendar defaults ---------------------------------------------------
  enrollment_details <- normalize_enrollment_details(enrollment_details)
  followup_details   <- normalize_followup_details(followup_details)
  trial_end_details  <- normalize_trial_end_details(trial_end_details)

  # ---- checks --------------------------------------------------------------
  # Note: update `check_makeData_args()` later to formally validate
  # `followup_details` and `trial_end_details`.
  chk <- check_makeData_args(
    correlation_matrix = correlation_matrix,
    SEED = SEED,
    sample_size_per_group = sample_size_per_group,
    endpoint_details = endpoint_details,
    enrollment_details = enrollment_details,
    followup_details = followup_details,
    trial_end_details = trial_end_details,
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

  # ---- design skeleton -----------------------------------------------------
  tt <- data.frame(
    trt = 0:(K - 1L),
    n_i = as.integer(n_by_arm)
  )

  # ---- helpers -------------------------------------------------------------
  expand_eff_K <- function(eff, K) {
    if (is.null(eff)) return(rep(0, K))
    if (!is.numeric(eff) || anyNA(eff)) {
      stop("`trt_effect` must be numeric, or NULL.")
    }
    if (length(eff) == 1L) return(c(0, rep(eff, K - 1L)))
    if (length(eff) == (K - 1L)) return(c(0, eff))
    stop("`trt_effect` must have length 1 or K-1, or be NULL.")
  }

  expand_param_K_including_control <- function(x, K, name = "param") {
    if (is.null(x)) stop("Internal error: missing ", name)
    if (!is.numeric(x) || anyNA(x)) {
      stop("`", name, "` must be numeric with no NA.")
    }
    if (length(x) == 1L) return(rep(x, K))
    if (length(x) == K) return(as.numeric(x))
    stop("`", name, "` must have length 1 or K. K = ", K, ".")
  }

  expand_active <- function(x, K, name = "arg") {
    if (is.null(x)) return(NULL)
    if (!is.numeric(x) || anyNA(x)) {
      stop("`", name, "` must be numeric with no NA.")
    }
    if (length(x) == 1L) return(rep(x, K - 1L))
    if (length(x) == (K - 1L)) return(as.numeric(x))
    stop("`", name, "` must have length 1 or K-1. K = ", K, ".")
  }

  # ---- build arm-specific marginal parameters ------------------------------
  for (j in seq_len(p)) {
    spec <- endpoint_details[[j]]
    typ  <- endpoint_types[j]

    if (typ == "continuous") {
      eff_vec <- expand_eff_K(spec$trt_effect %||% NULL, K)
      arm_eff <- eff_vec[tt$trt + 1L]
      tt[[paste0("mu", j)]] <- spec$baseline_mean + arm_eff

      sd_vec <- expand_param_K_including_control(spec$sd, K, name = "sd")
      if (any(sd_vec <= 0)) {
        stop("Continuous endpoint j = ", j, " `sd` must be > 0.")
      }
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

    } else {
      stop("Unknown endpoint type at j = ", j)
    }
  }

  if (!is.null(SEED)) set.seed(SEED)

  # ---- correlation handler -------------------------------------------------
  single_endpoint_mode <- is.null(correlation_matrix)

  if (!single_endpoint_mode) {
    L_latent_default <- make_pd(correlation_matrix)$L
  }

  cc <- calibration_control
  cc_n_mc      <- cc$n_mc           %||% 10000
  cc_tol       <- cc$tol            %||% 0.001
  cc_maxit     <- cc$maxit          %||% 100
  cc_rho_cap   <- cc$rho_cap        %||% 0.999
  cc_ensure_pd <- cc$ensure_pd      %||% TRUE
  cc_conv_type <- cc$conv_norm_type %||% "O"

  # ---- simulate per arm ----------------------------------------------------
  conditional_dist_sim <- function(i) {
    n_i <- tt$n_i[i]

    create_dist_function <- function(j, mu_value, sd_value = NULL) {
      spec <- endpoint_details[[j]]
      typ  <- endpoint_types[j]

      if (typ == "continuous") {
        local({
          m <- mu_value
          s <- sd_value
          function(u) stats::qnorm(u, mean = m, sd = s)
        })

      } else if (typ == "binary") {
        local({
          p <- mu_value
          function(u) stats::qbinom(u, size = 1, prob = p)
        })

      } else if (typ == "count") {
        local({
          mu <- mu_value
          size <- spec$size
          p0 <- spec$p_zero %||% 0
          function(u) qzinb_mixture(u, mu = mu, size = size, p0 = p0)
        })

      } else if (typ == "time-to-event") {
        local({
          r <- mu_value
          function(u) stats::qexp(u, rate = r)
        })

      } else {
        stop("Unknown endpoint type at j = ", j)
      }
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
      U <- matrix(stats::runif(n_i), ncol = 1)
      eps <- 1e-12
      U <- pmin(pmax(U, eps), 1 - eps)

    } else {
      if (isTRUE(target_correlation)) {
        latent_cor_cal <- calibrate_latent_cor_matrix(
          target_cor         = correlation_matrix,
          qfuns              = qfun_list,
          ensure_pd          = isTRUE(cc_ensure_pd),
          conv_norm_type     = cc_conv_type,
          return_diagnostics = FALSE,
          n_mc               = cc_n_mc,
          seed               = NULL,
          tol                = cc_tol,
          maxit              = cc_maxit,
          rho_cap            = cc_rho_cap
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

    if (!control_only) {
      sim_data$trt <- tt[i, "trt"]
    }

    sim_data
  }

  total_sim_data <- do.call(
    rbind,
    lapply(seq_len(nrow(tt)), conditional_dist_sim)
  )

  # ---- TTE independent censoring + indicators ------------------------------
  tte_idx <- which(endpoint_types == "time-to-event")

  if (length(tte_idx) > 0L) {
    censoring_rates <- vapply(
      tte_idx,
      function(j) endpoint_details[[j]]$censoring_rate %||% 0,
      numeric(1)
    )

    fatal_events <- vapply(
      tte_idx,
      function(j) isTRUE(endpoint_details[[j]]$fatal_event %||% FALSE),
      logical(1)
    )

    for (k in seq_along(tte_idx)) {
      j <- tte_idx[k]
      time_col <- paste0("V", j)

      ev_times <- total_sim_data[[time_col]]
      cr <- censoring_rates[k]

      cens_times <- if (cr <= 0) {
        rep(Inf, length(ev_times))
      } else {
        stats::rexp(length(ev_times), rate = cr)
      }

      is_censored <- ev_times > cens_times

      total_sim_data[[time_col]] <- pmin(ev_times, cens_times)
      total_sim_data[[paste0("Status_", k)]] <- as.integer(!is_censored)
    }

    # ---- semi-competing risks: non-fatal censoring can censor other TTEs ----
    if (non_fatal_censors_fatal) {
      fatal_k <- which(fatal_events)
      nonfatal_k <- setdiff(seq_along(tte_idx), fatal_k)

      for (k in nonfatal_k) {
        time_col_k <- paste0("V", tte_idx[k])
        cens_col_k <- paste0("Status_", k)

        is_cens <- total_sim_data[[cens_col_k]] == 0L
        if (!any(is_cens)) next

        t_cens <- total_sim_data[[time_col_k]][is_cens]

        for (kk in seq_along(tte_idx)) {
          if (kk == k) next

          time_col_kk <- paste0("V", tte_idx[kk])
          cens_col_kk <- paste0("Status_", kk)

          total_sim_data[[time_col_kk]][is_cens] <- pmin(
            total_sim_data[[time_col_kk]][is_cens],
            t_cens
          )

          total_sim_data[[cens_col_kk]][is_cens] <- 0L
        }
      }
    }

    # ---- fatal-event logic -------------------------------------------------
    if (length(tte_idx) > 1L && any(fatal_events)) {
      if (length(tte_idx) >= 2L && fatal_events[1] && fatal_events[2]) {
        t1 <- total_sim_data[[paste0("V", tte_idx[1])]]
        t2 <- total_sim_data[[paste0("V", tte_idx[2])]]
        c1 <- total_sim_data[["Status_1"]]
        c2 <- total_sim_data[["Status_2"]]

        cond1 <- (c1 == 1L) & (t1 < t2)
        total_sim_data[[paste0("V", tte_idx[2])]][cond1] <- t1[cond1]
        total_sim_data[["Status_2"]][cond1] <- 0L

        cond2 <- (c2 == 1L) & (t2 < t1)
        total_sim_data[[paste0("V", tte_idx[1])]][cond2] <- t2[cond2]
        total_sim_data[["Status_1"]][cond2] <- 0L

        cond3 <- (c1 == 0L) & (c2 == 0L)
        m <- pmin(t1[cond3], t2[cond3])

        total_sim_data[[paste0("V", tte_idx[1])]][cond3] <- m
        total_sim_data[[paste0("V", tte_idx[2])]][cond3] <- m
      }

      for (k in seq_along(tte_idx)) {
        if (!fatal_events[k]) next

        t_k <- total_sim_data[[paste0("V", tte_idx[k])]]

        if (k < length(tte_idx)) {
          for (kk in (k + 1L):length(tte_idx)) {
            time_col_kk <- paste0("V", tte_idx[kk])
            cens_col_kk <- paste0("Status_", kk)

            t_kk <- total_sim_data[[time_col_kk]]
            cond <- t_kk > t_k

            total_sim_data[[time_col_kk]][cond] <- t_k[cond]
            total_sim_data[[cens_col_kk]][cond] <- 0L
          }
        }
      }
    }
  }

  # ---- Enrollment, trial end, and administrative follow-up -----------------
  nSubs <- nrow(total_sim_data)

  enroll_time <- simulate_enrollment_times(
    n = nSubs,
    enrollment_details = enrollment_details,
    randomize_order = TRUE
  )

  add_enroll_column <-
    enrollment_details$enrollment_distribution != "none" ||
    trial_end_details$type != "none" ||
    !is.null(followup_details$max_followup)

  if (isTRUE(add_enroll_column)) {
    total_sim_data$enrollTime <- enroll_time
  }

  trial_calendar <- determine_trial_end_time(
    enroll_time = enroll_time,
    followup_details = followup_details,
    trial_end_details = trial_end_details,
    total_sim_data = total_sim_data,
    tte_idx = tte_idx
  )

  available_followup <- derive_available_followup(
    enroll_time = enroll_time,
    trial_end_time = trial_calendar$trial_end_time,
    followup_details = followup_details
  )

  add_followup_column <-
    trial_end_details$type != "none" ||
    !is.null(followup_details$max_followup)

  if (isTRUE(add_followup_column)) {
    total_sim_data$availableFollowup <- available_followup
  }

  if (length(tte_idx) > 0L && any(is.finite(available_followup))) {
    total_sim_data <- apply_followup_censoring(
      total_sim_data = total_sim_data,
      tte_idx = tte_idx,
      available_followup = available_followup
    )
  }

  # ---- rename endpoints ----------------------------------------------------
  cont_k <- 0L
  bin_k  <- 0L
  tte_k  <- 0L
  cnt_k  <- 0L

  new_names <- character(p)

  for (j in seq_len(p)) {
    if (endpoint_types[j] == "continuous") {
      cont_k <- cont_k + 1L
      new_names[j] <- paste0("Cont_", cont_k)
    }

    if (endpoint_types[j] == "binary") {
      bin_k <- bin_k + 1L
      new_names[j] <- paste0("Bin_", bin_k)
    }

    if (endpoint_types[j] == "count") {
      cnt_k <- cnt_k + 1L
      new_names[j] <- paste0("Int_", cnt_k)
    }

    if (endpoint_types[j] == "time-to-event") {
      tte_k <- tte_k + 1L
      new_names[j] <- paste0("TTE_", tte_k)
    }
  }

  names(total_sim_data)[
    match(paste0("V", seq_len(p)), names(total_sim_data))
  ] <- new_names

  # ---- metadata ------------------------------------------------------------
  meta <- list(
    correlation_matrix   = correlation_matrix,
    target_correlation   = isTRUE(target_correlation) && !single_endpoint_mode,
    endpoint_details     = endpoint_details,
    endpoint_types       = endpoint_types,
    endpoint_names       = new_names,
    n_arms               = K,
    n_by_arm             = as.integer(n_by_arm),
    control_only         = control_only,
    single_endpoint_mode = single_endpoint_mode,
    enrollment_details   = enrollment_details,
    followup_details     = followup_details,
    trial_end_details    = trial_end_details,
    trial_calendar       = trial_calendar
  )

  new_makeDataSim(total_sim_data, meta)
}
