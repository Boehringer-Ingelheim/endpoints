# util_funcs.R ---------------------------------------------------------------

# helper convention used throughout: if a is not NULL return a, else return b
`%||%` <- function(a, b) if (!is.null(a)) a else b


# -----------------------------------------------------------------------------
# Basic math helpers
# -----------------------------------------------------------------------------

logit <- function(p) {
  if (anyNA(p)) stop("logit(): p contains NA.")
  if (any(p <= 0 | p >= 1)) stop("logit(): p must be in (0, 1).")
  stats::qlogis(p)
}

inv_logit <- function(x) stats::plogis(x)


# -----------------------------------------------------------------------------
# Rate helper
# -----------------------------------------------------------------------------

#' Solve exponential rate parameters from target observed event probabilities
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

  if (mode == "simple") {
    if (!is.numeric(event_rate) || length(event_rate) != 1L ||
        is.na(event_rate) || event_rate <= 0) {
      stop("For mode = 'simple', `event_rate` must be a single positive number.")
    }

    return(event_rate / target_prob - event_rate)
  }

  if (mode == "admin") {
    if (!is.numeric(admin_time) || length(admin_time) != 1L ||
        is.na(admin_time) || admin_time <= 0) {
      stop("For mode = 'admin', `admin_time` must be a single positive number.")
    }

    return(-log(1 - target_prob) / admin_time)
  }

  if (mode == "semi-competing") {
    vals <- c(fatal_event_rate, fatal_censor_rate, nonfatal_event_rate)

    if (!is.numeric(vals) || anyNA(vals) || any(vals < 0)) {
      stop(
        "For mode = 'semi-competing', `fatal_event_rate`, ",
        "`fatal_censor_rate`, and `nonfatal_event_rate` must be numeric ",
        "and >= 0."
      )
    }

    if (length(fatal_event_rate) != 1L ||
        length(fatal_censor_rate) != 1L ||
        length(nonfatal_event_rate) != 1L) {
      stop(
        "For mode = 'semi-competing', `fatal_event_rate`, ",
        "`fatal_censor_rate`, and `nonfatal_event_rate` must each be length 1."
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


# -----------------------------------------------------------------------------
# Correlation matrix helper
# -----------------------------------------------------------------------------

#' Construct a correlation matrix from endpoint-index triplets
#'
#' @export
corr_make <- function(num_endpoints, values = NULL) {
  if (!is.numeric(num_endpoints) || length(num_endpoints) != 1L ||
      is.na(num_endpoints)) {
    stop("`num_endpoints` must be a single positive integer.")
  }

  num_endpoints <- as.integer(num_endpoints)

  if (num_endpoints < 1L) {
    stop("`num_endpoints` must be >= 1.")
  }

  R <- diag(num_endpoints)

  if (is.null(values)) {
    return(R)
  }

  vals <- matrix(as.matrix(values), ncol = 3)

  if (nrow(vals) == 0L) {
    return(R)
  }

  if (anyNA(vals)) {
    stop("`values` must not contain NA.")
  }

  ij <- vals[, 1:2, drop = FALSE]
  rho <- vals[, 3]

  if (any(ij <= 0)) {
    stop("Indices in `values[, 1:2]` must be positive.")
  }

  if (any(ij > num_endpoints)) {
    stop("Indices in `values[, 1:2]` exceed `num_endpoints`.")
  }

  if (any(rho < -1 | rho > 1)) {
    stop("All correlations in `values[, 3]` must be in [-1, 1].")
  }

  i <- pmin.int(ij[, 1], ij[, 2])
  j <- pmax.int(ij[, 1], ij[, 2])

  R[cbind(i, j)] <- rho
  R[cbind(j, i)] <- rho

  diag(R) <- 1

  R
}


# -----------------------------------------------------------------------------
# Endpoint type helper
# -----------------------------------------------------------------------------

normalize_endpoint_type <- function(x) {
  x <- tolower(x)

  if (x %in% c("normal", "continuous", "gaussian")) {
    return("continuous")
  }

  if (x %in% c("binary", "bin")) {
    return("binary")
  }

  if (x %in% c("count", "nb", "negbinom", "negative binomial", "zinb")) {
    return("count")
  }

  if (x %in% c("tte", "time-to-event", "time_to_event", "survival")) {
    return("time-to-event")
  }

  stop("Unknown endpoint_type: ", x)
}


# -----------------------------------------------------------------------------
# Positive-definite correlation helpers
# -----------------------------------------------------------------------------

make_pd <- function(R, eps = 1e-10) {
  out <- tryCatch(chol(R), error = function(e) NULL)

  if (!is.null(out)) {
    return(list(R = R, L = out))
  }

  ee <- eigen(R, symmetric = TRUE)
  vals2 <- pmax(ee$values, eps)

  R2 <- ee$vectors %*% diag(vals2, nrow = length(vals2)) %*% t(ee$vectors)

  d <- sqrt(diag(R2))
  R2 <- sweep(sweep(R2, 1, d, "/"), 2, d, "/")

  L2 <- chol(R2)

  list(R = R2, L = L2)
}


project_to_pd_cor <- function(R, conv_norm_type = "F") {
  if (requireNamespace("Matrix", quietly = TRUE)) {
    pmt <- Matrix::nearPD(
      R,
      corr = TRUE,
      keepDiag = TRUE,
      conv.norm.type = conv_norm_type,
      trace = FALSE
    )

    return(as.matrix(pmt$mat))
  }

  make_pd(R)$R
}


# -----------------------------------------------------------------------------
# Gaussian copula
# -----------------------------------------------------------------------------

draw_gaussian_copula_u <- function(n, L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  p <- ncol(L)

  Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p) %*% L
  U <- stats::pnorm(Z)

  eps <- 1e-12
  pmin(pmax(U, eps), 1 - eps)
}


# -----------------------------------------------------------------------------
# Marginal quantile helpers
# -----------------------------------------------------------------------------

qzinb_mixture <- function(u, mu, size, p0) {
  out <- integer(length(u))

  is0 <- u < p0
  out[is0] <- 0L

  if (any(!is0)) {
    u2 <- (u[!is0] - p0) / (1 - p0)
    out[!is0] <- stats::qnbinom(u2, size = size, mu = mu)
  }

  out
}


# -----------------------------------------------------------------------------
# Enrollment, follow-up, and trial-calendar helpers
# -----------------------------------------------------------------------------

normalize_enrollment_details <- function(enrollment_details = list()) {
  if (!is.list(enrollment_details)) {
    stop("`enrollment_details` must be a list.")
  }

  if ("administrative_censoring" %in% names(enrollment_details)) {
    warning(
      "`enrollment_details$administrative_censoring` is deprecated and ignored. ",
      "Use `followup_details` and `trial_end_details` instead."
    )

    enrollment_details$administrative_censoring <- NULL
  }

  enrollment_details <- utils::modifyList(
    list(
      enrollment_distribution        = "none",
      enrollment_exponential_rate    = NULL,
      piecewise_enrollment_cutpoints = NULL,
      piecewise_enrollment_rates     = NULL
    ),
    enrollment_details
  )

  enrollment_details$enrollment_distribution <- match.arg(
    enrollment_details$enrollment_distribution,
    c("none", "exponential", "piecewise")
  )

  enrollment_details
}


normalize_followup_details <- function(followup_details = list()) {
  if (!is.list(followup_details)) {
    stop("`followup_details` must be a list.")
  }

  followup_details <- utils::modifyList(
    list(
      min_followup = NULL,
      max_followup = NULL
    ),
    followup_details
  )

  if (!is.null(followup_details$min_followup)) {
    if (!is.numeric(followup_details$min_followup) ||
        length(followup_details$min_followup) != 1L ||
        is.na(followup_details$min_followup) ||
        followup_details$min_followup < 0) {
      stop("`followup_details$min_followup` must be a non-negative numeric scalar.")
    }
  }

  if (!is.null(followup_details$max_followup)) {
    if (!is.numeric(followup_details$max_followup) ||
        length(followup_details$max_followup) != 1L ||
        is.na(followup_details$max_followup) ||
        followup_details$max_followup <= 0) {
      stop("`followup_details$max_followup` must be a positive numeric scalar.")
    }
  }

  if (!is.null(followup_details$min_followup) &&
      !is.null(followup_details$max_followup) &&
      followup_details$min_followup > followup_details$max_followup) {
    warning(
      "`followup_details$min_followup` is greater than ",
      "`followup_details$max_followup`. Check that this is intentional and ",
      "that time units are consistent."
    )
  }

  followup_details
}


normalize_trial_end_details <- function(trial_end_details = list()) {
  if (!is.list(trial_end_details)) {
    stop("`trial_end_details` must be a list.")
  }

  trial_end_details <- utils::modifyList(
    list(
      type                 = "none",
      trial_end_time       = NULL,
      event_endpoint       = NULL,
      target_events        = NULL,
      require_min_followup = FALSE,
      max_trial_duration   = NULL,
      target_not_reached   = "error"
    ),
    trial_end_details
  )

  trial_end_details$type <- match.arg(
    trial_end_details$type,
    c("none", "fixed_calendar", "last_patient_min_followup", "event_driven")
  )

  trial_end_details$target_not_reached <- match.arg(
    trial_end_details$target_not_reached,
    c("error", "max_trial_duration", "last_patient_min_followup")
  )

  if (!is.logical(trial_end_details$require_min_followup) ||
      length(trial_end_details$require_min_followup) != 1L ||
      is.na(trial_end_details$require_min_followup)) {
    stop("`trial_end_details$require_min_followup` must be TRUE or FALSE.")
  }

  trial_end_details
}


make_poisson_enrollment <- function(n, rate) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
    stop("`n` must be a positive integer.")
  }

  if (!is.numeric(rate) || length(rate) != 1L || is.na(rate) || rate <= 0) {
    stop("`rate` must be a positive numeric scalar.")
  }

  gaps <- stats::rexp(n, rate = rate)
  enroll_time <- cumsum(gaps)

  # Start the trial clock at first patient randomized.
  enroll_time - min(enroll_time)
}


make_piecewise_poisson_enrollment <- function(n, cutpoints, rates) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
    stop("`n` must be a positive integer.")
  }

  n <- as.integer(n)

  if (!is.numeric(cutpoints) || length(cutpoints) < 2L || anyNA(cutpoints)) {
    stop("`piecewise_enrollment_cutpoints` must be numeric with at least two values.")
  }

  if (is.unsorted(cutpoints, strictly = TRUE)) {
    stop("`piecewise_enrollment_cutpoints` must be strictly increasing.")
  }

  if (cutpoints[1] != 0) {
    stop("`piecewise_enrollment_cutpoints` should start at 0.")
  }

  if (!is.numeric(rates) || anyNA(rates) || any(rates < 0)) {
    stop("`piecewise_enrollment_rates` must be non-negative numeric values.")
  }

  if (length(rates) != length(cutpoints) - 1L) {
    stop("`piecewise_enrollment_rates` must have length `length(cutpoints) - 1`.")
  }

  if (all(rates == 0)) {
    stop("At least one piecewise enrollment rate must be > 0.")
  }

  enroll_time <- numeric(n)
  n_enrolled <- 0L

  t <- cutpoints[1]
  k <- 1L

  while (n_enrolled < n) {
    # If we have moved beyond the user-specified intervals, continue with the
    # final interval rate.
    if (k > length(rates)) {
      final_rate <- tail(rates, 1L)

      if (final_rate <= 0) {
        stop(
          "The specified piecewise enrollment process did not generate enough ",
          "subjects, and the final enrollment rate is 0. Increase the final rate ",
          "or extend the cutpoints."
        )
      }

      gap <- stats::rexp(1L, rate = final_rate)
      t <- t + gap

      n_enrolled <- n_enrolled + 1L
      enroll_time[n_enrolled] <- t

      next
    }

    interval_end <- cutpoints[k + 1L]
    rate_k <- rates[k]

    # If the current interval has no accrual, jump to the next interval.
    if (rate_k <= 0) {
      t <- max(t, interval_end)
      k <- k + 1L
      next
    }

    # Draw waiting time under the current interval's rate.
    gap <- stats::rexp(1L, rate = rate_k)
    proposed_time <- t + gap

    if (proposed_time < interval_end) {
      # Arrival occurs inside the current interval.
      t <- proposed_time
      n_enrolled <- n_enrolled + 1L
      enroll_time[n_enrolled] <- t
    } else {
      # No arrival before the interval boundary. Move to the next interval.
      #
      # Memorylessness makes this valid: once we reach the new interval,
      # the waiting time is redrawn using the new rate.
      t <- interval_end
      k <- k + 1L
    }
  }

  # Start the trial clock at first patient randomized.
  enroll_time - min(enroll_time)
}


simulate_enrollment_times <- function(n,
                                      enrollment_details,
                                      randomize_order = TRUE) {
  dist <- enrollment_details$enrollment_distribution

  if (dist == "none") {
    enroll_time <- rep(0, n)

  } else if (dist == "exponential") {
    rate <- enrollment_details$enrollment_exponential_rate

    if (is.null(rate) || !is.numeric(rate) || length(rate) != 1L ||
        is.na(rate) || rate <= 0) {
      stop(
        "For `enrollment_distribution = 'exponential'`, provide ",
        "`enrollment_exponential_rate > 0`."
      )
    }

    enroll_time <- make_poisson_enrollment(n = n, rate = rate)

  } else if (dist == "piecewise") {
    enroll_time <- make_piecewise_poisson_enrollment(
      n         = n,
      cutpoints = enrollment_details$piecewise_enrollment_cutpoints,
      rates     = enrollment_details$piecewise_enrollment_rates
    )

  } else {
    stop("Unknown enrollment distribution: ", dist)
  }

  # makeData() creates rows by arm, so randomize enrollment assignment to avoid
  # confounding arm with calendar enrollment order.
  if (isTRUE(randomize_order)) {
    enroll_time <- sample(enroll_time, size = n, replace = FALSE)
  }

  enroll_time
}


resolve_event_endpoint <- function(event_endpoint, tte_idx) {
  if (length(tte_idx) == 0L) {
    stop("Event-driven trial end requires at least one TTE endpoint.")
  }

  if (is.null(event_endpoint)) {
    stop("For event-driven trial end, provide `trial_end_details$event_endpoint`.")
  }

  if (is.numeric(event_endpoint) && length(event_endpoint) == 1L) {
    k <- as.integer(event_endpoint)

    if (is.na(k) || k < 1L || k > length(tte_idx)) {
      stop(
        "Numeric `event_endpoint` is interpreted as TTE ordinal and must be ",
        "between 1 and the number of TTE endpoints."
      )
    }

    return(k)
  }

  if (is.character(event_endpoint) && length(event_endpoint) == 1L) {
    x <- event_endpoint

    if (grepl("^TTE_[0-9]+$", x)) {
      k <- as.integer(sub("^TTE_", "", x))
      if (k >= 1L && k <= length(tte_idx)) return(k)
    }

    if (grepl("^Status_[0-9]+$", x)) {
      k <- as.integer(sub("^Status_", "", x))
      if (k >= 1L && k <= length(tte_idx)) return(k)
    }

    if (grepl("^V[0-9]+$", x)) {
      j <- as.integer(sub("^V", "", x))
      k <- match(j, tte_idx)
      if (!is.na(k)) return(k)
    }
  }

  stop(
    "`event_endpoint` must be a TTE ordinal such as 1, or a character value ",
    "such as 'TTE_1', 'Status_1', or 'Vj'."
  )
}


determine_trial_end_time <- function(enroll_time,
                                     followup_details,
                                     trial_end_details,
                                     total_sim_data = NULL,
                                     tte_idx = integer(0),
                                     tol = 1e-10) {
  type <- trial_end_details$type

  min_followup <- followup_details$min_followup
  max_followup <- followup_details$max_followup

  if (!is.numeric(enroll_time) || anyNA(enroll_time)) {
    stop("`enroll_time` must be numeric with no NA.")
  }

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol < 0) {
    stop("`tol` must be a non-negative numeric scalar.")
  }

  lpi_min_followup_time <- NULL

  if (!is.null(min_followup)) {
    lpi_min_followup_time <- max(enroll_time) + min_followup
  }

  out <- list(
    trial_end_time       = Inf,
    trial_end_reason     = "none",
    event_target_reached = NA,
    n_events_at_end      = NA_integer_
  )

  if (type == "none") {
    return(out)
  }

  if (type == "fixed_calendar") {
    trial_end_time <- trial_end_details$trial_end_time

    if (is.null(trial_end_time) || !is.numeric(trial_end_time) ||
        length(trial_end_time) != 1L || is.na(trial_end_time) ||
        trial_end_time < 0) {
      stop(
        "For `trial_end_details$type = 'fixed_calendar'`, provide ",
        "`trial_end_details$trial_end_time >= 0`."
      )
    }

    out$trial_end_time <- trial_end_time
    out$trial_end_reason <- "fixed_calendar"

    return(out)
  }

  if (type == "last_patient_min_followup") {
    if (is.null(min_followup)) {
      stop(
        "`followup_details$min_followup` is required when ",
        "`trial_end_details$type = 'last_patient_min_followup'`."
      )
    }

    out$trial_end_time <- lpi_min_followup_time
    out$trial_end_reason <- "last_patient_min_followup"

    return(out)
  }

  if (type == "event_driven") {
    if (is.null(total_sim_data)) {
      stop("`total_sim_data` is required for event-driven trial end.")
    }

    target_events <- trial_end_details$target_events

    if (is.null(target_events) || !is.numeric(target_events) ||
        length(target_events) != 1L || is.na(target_events) ||
        target_events < 1) {
      stop(
        "For `trial_end_details$type = 'event_driven'`, provide ",
        "`trial_end_details$target_events` as a positive integer."
      )
    }

    target_events <- as.integer(target_events)

    k <- resolve_event_endpoint(
      event_endpoint = trial_end_details$event_endpoint,
      tte_idx = tte_idx
    )

    time_col <- paste0("V", tte_idx[k])
    status_col <- paste0("Status_", k)

    if (!all(c(time_col, status_col) %in% names(total_sim_data))) {
      stop(
        "Could not find required TTE columns for event-driven trial end: ",
        time_col, " and ", status_col, "."
      )
    }

    tte_time <- total_sim_data[[time_col]]
    tte_status <- total_sim_data[[status_col]]

    if (!is.numeric(tte_time) || anyNA(tte_time)) {
      stop("Event-driven TTE time column `", time_col, "` must be numeric with no NA.")
    }

    if (!is.numeric(tte_status) && !is.integer(tte_status)) {
      stop("Event-driven status column `", status_col, "` must be numeric/integer.")
    }

    if (length(tte_time) != length(enroll_time)) {
      stop("Length of `enroll_time` must match the number of rows in `total_sim_data`.")
    }

    eligible <- tte_status == 1L

    # If max_followup is a subject-level cap, events after max_followup would
    # not be observable and should not drive trial end.
    #
    # Use a tolerance so events exactly at max_followup are still considered
    # observable despite floating-point noise.
    if (!is.null(max_followup)) {
      eligible <- eligible & tte_time <= max_followup + tol
    }

    calendar_event_times <- enroll_time[eligible] + tte_time[eligible]
    calendar_event_times <- sort(calendar_event_times)

    target_reached <- length(calendar_event_times) >= target_events

    if (target_reached) {
      trial_end_time <- calendar_event_times[target_events]
      trial_end_reason <- "event_target"

    } else {
      if (!is.null(trial_end_details$max_trial_duration) &&
          trial_end_details$target_not_reached == "max_trial_duration") {
        trial_end_time <- trial_end_details$max_trial_duration
        trial_end_reason <- "max_trial_duration_target_not_reached"

        warning(
          "Event target was not reached. Using ",
          "`trial_end_details$max_trial_duration` as trial end."
        )

      } else if (!is.null(lpi_min_followup_time) &&
                 trial_end_details$target_not_reached == "last_patient_min_followup") {
        trial_end_time <- lpi_min_followup_time
        trial_end_reason <- "last_patient_min_followup_target_not_reached"

        warning(
          "Event target was not reached. Using last-patient-min-follow-up ",
          "time as trial end."
        )

      } else {
        stop(
          "Event target was not reached. Consider increasing sample size, ",
          "event rate, follow-up, or specifying `max_trial_duration` with ",
          "`target_not_reached = 'max_trial_duration'`."
        )
      }
    }

    if (isTRUE(trial_end_details$require_min_followup)) {
      if (is.null(lpi_min_followup_time)) {
        stop(
          "`trial_end_details$require_min_followup = TRUE` requires ",
          "`followup_details$min_followup`."
        )
      }

      trial_end_time <- max(trial_end_time, lpi_min_followup_time)
      trial_end_reason <- paste0(trial_end_reason, "_plus_min_followup")
    }

    if (!is.null(trial_end_details$max_trial_duration)) {
      max_trial_duration <- trial_end_details$max_trial_duration

      if (!is.numeric(max_trial_duration) ||
          length(max_trial_duration) != 1L ||
          is.na(max_trial_duration) ||
          max_trial_duration <= 0) {
        stop("`trial_end_details$max_trial_duration` must be a positive numeric scalar.")
      }

      if (trial_end_time > max_trial_duration) {
        trial_end_time <- max_trial_duration
        trial_end_reason <- "max_trial_duration"
      }
    }

    out$trial_end_time <- trial_end_time
    out$trial_end_reason <- trial_end_reason
    out$event_target_reached <- target_reached

    # Count all events occurring at or before the final trial end time.
    # The tolerance prevents the event defining the trial end from being missed
    # due to tiny floating-point differences.
    out$n_events_at_end <- if (length(calendar_event_times) == 0L) {
      0L
    } else {
      sum(calendar_event_times <= trial_end_time + tol)
    }

    return(out)
  }

  stop("Unknown trial end type: ", type)
}


derive_available_followup <- function(enroll_time,
                                      trial_end_time,
                                      followup_details) {
  if (is.infinite(trial_end_time)) {
    available_followup <- rep(Inf, length(enroll_time))
  } else {
    available_followup <- pmax(0, trial_end_time - enroll_time)
  }

  if (!is.null(followup_details$max_followup)) {
    available_followup <- pmin(
      available_followup,
      followup_details$max_followup
    )
  }

  available_followup
}



apply_followup_censoring <- function(total_sim_data,
                                     tte_idx,
                                     available_followup,
                                     tol = 1e-10) {
  if (length(tte_idx) == 0L) {
    return(total_sim_data)
  }

  if (!is.numeric(available_followup) || anyNA(available_followup)) {
    stop("`available_followup` must be numeric with no NA.")
  }

  for (k in seq_along(tte_idx)) {
    time_col <- paste0("V", tte_idx[k])
    status_col <- paste0("Status_", k)

    t <- total_sim_data[[time_col]]

    # Administrative censoring should only happen when the current observed
    # time is meaningfully beyond available follow-up.
    #
    # Events exactly at the cutoff should remain observed events.
    is_admin_cens <- t > available_followup + tol

    total_sim_data[[time_col]] <- pmin(t, available_followup)
    total_sim_data[[status_col]][is_admin_cens] <- 0L
  }

  total_sim_data
}



# -----------------------------------------------------------------------------
# Correlation calibration
# -----------------------------------------------------------------------------

calibrate_latent_rho_pair <- function(target_r,
                                      qfun1,
                                      qfun2,
                                      n_mc = 10000,
                                      seed = NULL,
                                      tol = 0.001,
                                      maxit = 100,
                                      rho_cap = 0.999) {
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

  z1_base <- stats::rnorm(n_mc)
  z2_base <- stats::rnorm(n_mc)

  observed_corr_minus_target <- function(rho) {
    rho <- max(min(rho, rho_cap), -rho_cap)

    z1 <- z1_base
    z2 <- rho * z1_base + sqrt(pmax(0, 1 - rho^2)) * z2_base

    u1 <- stats::pnorm(z1)
    u2 <- stats::pnorm(z2)

    x1 <- qfun1(u1)
    x2 <- qfun2(u2)

    if (stats::sd(x1) == 0 || stats::sd(x2) == 0) {
      return(NA_real_)
    }

    stats::cor(x1, x2) - target_r
  }

  lo <- -rho_cap
  hi <-  rho_cap

  f_lo <- observed_corr_minus_target(lo)
  f_hi <- observed_corr_minus_target(hi)

  if (is.na(f_lo) || is.na(f_hi)) {
    return(list(root = 0, note = "degenerate margin"))
  }

  if (f_lo * f_hi > 0) {
    root_pick <- if (abs(f_lo) < abs(f_hi)) lo else hi
    return(list(root = root_pick, note = "target infeasible; picked closest endpoint"))
  }

  stats::uniroot(
    f = observed_corr_minus_target,
    interval = c(lo, hi),
    tol = tol,
    maxiter = maxit
  )
}


calibrate_latent_cor_matrix <- function(target_cor,
                                        qfuns,
                                        ensure_pd = TRUE,
                                        conv_norm_type = "F",
                                        return_diagnostics = FALSE,
                                        n_mc = 10000,
                                        seed = NULL,
                                        tol = 0.001,
                                        maxit = 100,
                                        rho_cap = 0.999) {
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

  if (p >= 2L) {
    for (r in seq_len(p - 1L)) {
      conv_diag[[r]] <- vector("list", p)

      for (c in (r + 1L):p) {
        seed_rc <- if (is.null(seed)) NULL else seed + 1000L * r + c

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
  cat("  n = ", nrow(dat), "\n", sep = "")
  cat("  n_arms = ", K, "\n", sep = "")
  cat("  endpoints = ", length(meta$endpoint_details %||% list()), "\n", sep = "")
  cat("  target_correlation = ", isTRUE(meta$target_correlation), "\n", sep = "")
  cat("  single_endpoint_mode = ", isTRUE(meta$single_endpoint_mode), "\n", sep = "")
  cat("  control_only = ", isTRUE(meta$control_only), "\n", sep = "")

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

  if (!is.null(meta$trial_calendar)) {
    tc <- meta$trial_calendar

    cat("  trial_end_reason = ", tc$trial_end_reason %||% "none", "\n", sep = "")

    if (!is.null(tc$trial_end_time)) {
      cat("  trial_end_time = ", tc$trial_end_time, "\n", sep = "")
    }

    if (!is.na(tc$event_target_reached %||% NA)) {
      cat("  event_target_reached = ", tc$event_target_reached, "\n", sep = "")
    }

    if (!is.na(tc$n_events_at_end %||% NA)) {
      cat("  n_events_at_end = ", tc$n_events_at_end, "\n", sep = "")
    }
  }

  show_rows <- 6L
  show_cols <- 10L

  disp <- dat
  n <- nrow(disp)
  p <- ncol(disp)

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
    Xa <- if (has_trt) {
      dat[dat$trt == a, endpoint_names, drop = FALSE]
    } else {
      dat[, endpoint_names, drop = FALSE]
    }

    cor_by_arm[[paste0("arm_", a)]] <- round(
      suppressWarnings(stats::cor(Xa)),
      3
    )
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

  if (length(idx_cont) > 0L) {
    cont_tbl <- do.call(rbind, lapply(idx_cont, function(j) {
      col  <- endpoint_names[j]
      spec <- endpoint_details[[j]]

      fit <- if (has_trt && K > 1L) {
        stats::lm(dat[[col]] ~ trtF)
      } else {
        stats::lm(dat[[col]] ~ 1)
      }

      b0 <- unname(stats::coef(fit)[1])

      sd_in <- spec$sd

      sd_in_vec <- if (length(sd_in) == 1L) {
        rep(sd_in, K)
      } else if (length(sd_in) == K) {
        as.numeric(sd_in)
      } else {
        rep(NA_real_, K)
      }

      est_eff <- rep(0, K)
      names(est_eff) <- 0:(K - 1L)

      if (has_trt && K > 1L) {
        cf <- stats::coef(fit)

        for (a in 1:(K - 1L)) {
          nm <- paste0("trtF", a)
          if (nm %in% names(cf)) est_eff[a + 1L] <- unname(cf[nm])
        }
      }

      inp_eff <- expand_input_effect(spec$trt_effect %||% NULL, K)

      est_sd_by_arm <- vapply(0:(K - 1L), function(a) {
        idx <- arm_subset(a)
        if (sum(idx) <= 1L) return(NA_real_)
        stats::sd(stats::residuals(fit)[idx])
      }, numeric(1))

      do.call(rbind, lapply(0:(K - 1L), function(a) {
        data.frame(
          endpoint = col,
          arm = a,
          input_baseline_mean = spec$baseline_mean,
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

  if (length(idx_bin) > 0L) {
    bin_tbl <- do.call(rbind, lapply(idx_bin, function(j) {
      col  <- endpoint_names[j]
      spec <- endpoint_details[[j]]

      fit <- if (has_trt && K > 1L) {
        stats::glm(dat[[col]] ~ trtF, family = stats::binomial())
      } else {
        stats::glm(dat[[col]] ~ 1, family = stats::binomial())
      }

      cf <- stats::coef(fit)
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

  if (length(idx_cnt) > 0L) {
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

      cf <- stats::coef(fit)
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
          input_p_zero        = spec$p_zero %||% 0,
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

  if (length(idx_tte) > 0L) {
    if (!requireNamespace("survival", quietly = TRUE)) {
      stop("survival is required for TTE endpoint summary via survival::coxph().")
    }

    tte_tbl <- do.call(rbind, lapply(seq_along(idx_tte), function(k) {
      j    <- idx_tte[k]
      col  <- endpoint_names[j]
      spec <- endpoint_details[[j]]

      censor_col <- paste0("Status_", k)

      if (!censor_col %in% names(dat)) {
        stop("Missing censoring column ", censor_col, " for ", col, ".")
      }

      time <- dat[[col]]
      status <- dat[[censor_col]]

      fit <- tryCatch(
        {
          if (has_trt && K > 1L) {
            survival::coxph(survival::Surv(time, status) ~ trtF)
          } else {
            survival::coxph(survival::Surv(time, status) ~ 1)
          }
        },
        error = function(e) NULL
      )

      est_eff <- rep(NA_real_, K)
      names(est_eff) <- 0:(K - 1L)
      est_eff[1L] <- 0

      if (!is.null(fit) && has_trt && K > 1L) {
        cf <- stats::coef(fit)

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
        sa <- status[idx]
        ta <- time[idx]

        denom <- sum(ta)

        if (denom <= 0) {
          return(NA_real_)
        }

        sum(sa) / denom
      }, numeric(1))

      do.call(rbind, lapply(0:(K - 1L), function(a) {
        data.frame(
          endpoint = col,
          arm = a,
          censor_col = censor_col,
          input_baseline_rate = spec$baseline_rate,
          input_trt_logHR     = inp_eff[a + 1L],
          input_trt_HR        = exp(inp_eff[a + 1L]),
          est_trt_logHR       = est_eff[a + 1L],
          est_trt_HR          = exp(est_eff[a + 1L]),
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
    n_arms = K,
    enrollment_details = meta$enrollment_details %||% NULL,
    followup_details   = meta$followup_details %||% NULL,
    trial_end_details  = meta$trial_end_details %||% NULL,
    trial_calendar     = meta$trial_calendar %||% NULL
  )

  class(out) <- "summary.makeDataSim"

  out
}


#' @exportS3Method
print.summary.makeDataSim <- function(x, ...) {
  cat("<summary.makeDataSim>\n\n")
  cat("n_arms:", x$n_arms, "\n\n")

  if (!is.null(x$trial_calendar)) {
    cat("Trial calendar:\n")
    print(x$trial_calendar)
    cat("\n")
  }

  cat("Target correlation:\n")
  print(x$target_correlation)

  cat("\nEstimated correlation (by arm):\n")

  for (nm in names(x$estimated_correlation_by_arm)) {
    cat("\n", nm, ":\n", sep = "")
    print(x$estimated_correlation_by_arm[[nm]])
  }

  if (!is.null(x$continuous)) {
    cat("\nContinuous endpoints (endpoint x arm):\n")
    print(x$continuous)
  }

  if (!is.null(x$binary)) {
    cat("\nBinary endpoints (endpoint x arm):\n")
    print(x$binary)
  }

  if (!is.null(x$count)) {
    cat("\nCount endpoints (endpoint x arm):\n")
    print(x$count)
  }

  if (!is.null(x$tte)) {
    cat("\nTTE endpoints (endpoint x arm):\n")
    print(x$tte)
  }

  invisible(x)
}


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
    if (!is.numeric(arm) || length(arm) != 1L || is.na(arm)) {
      stop("`arm` must be a single number.")
    }

    arm <- as.integer(arm)

    if (arm < 0 || arm > (K - 1L)) {
      stop("`arm` must be in {0,1,...,", K - 1L, "}.")
    }

    dat_sub <- dat[dat$trt == arm, meta$endpoint_names, drop = FALSE]
    arm_lbl <- arm
  }

  if (!is.null(names)) {
    if (!is.character(names)) {
      stop("`names` must be a character vector.")
    }

    if (length(names) != ncol(dat_sub)) {
      stop(
        "`names` must have length equal to the number of endpoints (",
        ncol(dat_sub),
        ")."
      )
    }

    colnames(dat_sub) <- names
  }

  GGally::ggpairs(
    dat_sub,
    upper = list(continuous = "points"),
    lower = list(continuous = "cor"),
    progress = FALSE,
    ...
  ) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0("Arm ", arm_lbl))
}


# -----------------------------------------------------------------------------
# Argument checks for makeData()
# -----------------------------------------------------------------------------

check_makeData_args <- function(correlation_matrix,
                                SEED,
                                sample_size_per_group,
                                endpoint_details,
                                enrollment_details,
                                followup_details,
                                trial_end_details,
                                non_fatal_censors_fatal,
                                target_correlation,
                                calibration_control,
                                arm_mode = c("auto", "full", "control")) {

  arm_mode <- match.arg(arm_mode)

  # ---- endpoint_details ----------------------------------------------------
  if (!is.list(endpoint_details) || length(endpoint_details) < 1L) {
    stop("Error: `endpoint_details` must be a non-empty list of endpoint spec lists.")
  }

  p <- length(endpoint_details)

  endpoint_types <- vapply(endpoint_details, function(e) {
    if (is.null(e$endpoint_type)) {
      stop("Error: Each endpoint spec must include `endpoint_type`.")
    }

    normalize_endpoint_type(e$endpoint_type)
  }, character(1))

  # ---- SEED ----------------------------------------------------------------
  if (!is.null(SEED)) {
    if (!is.numeric(SEED) || length(SEED) != 1L || is.na(SEED)) {
      stop("Error: `SEED` must be NULL or a single numeric value.")
    }
  }

  # ---- correlation_matrix --------------------------------------------------
  single_endpoint_mode <- is.null(correlation_matrix)

  if (single_endpoint_mode) {
    if (p != 1L) {
      stop(
        "Error: If `correlation_matrix = NULL`, you must supply exactly one ",
        "endpoint."
      )
    }
  } else {
    if (!is.matrix(correlation_matrix) ||
        nrow(correlation_matrix) != p ||
        ncol(correlation_matrix) != p) {
      stop(
        "Error: Dimensions of `correlation_matrix` must match the number of ",
        "endpoints: ", p, " x ", p, "."
      )
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

  # ---- infer number of arms ------------------------------------------------
  active_len <- function(x, nm) {
    if (is.null(x)) return(0L)

    if (!is.numeric(x) || anyNA(x)) {
      stop("Error: `", nm, "` must be numeric with no NA.")
    }

    as.integer(length(x))
  }

  lens_effect <- vapply(
    endpoint_details,
    function(e) active_len(e$trt_effect %||% NULL, "trt_effect"),
    integer(1)
  )

  lens_bprob <- vapply(
    endpoint_details,
    function(e) active_len(e$trt_prob %||% NULL, "trt_prob"),
    integer(1)
  )

  lens_cmean <- vapply(
    endpoint_details,
    function(e) active_len(e$trt_count %||% NULL, "trt_count"),
    integer(1)
  )

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
        stop(
          "Error: endpoint ", j, " `", nm, "` must have length 1 or K-1 ",
          "(K = ", K, ")."
        )
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
      stop(
        "Error: control-only mode prohibits specifying `trt_effect`, ",
        "`trt_prob`, or `trt_count`."
      )
    }
  }

  # ---- sample size ---------------------------------------------------------
  if (!is.numeric(sample_size_per_group) || anyNA(sample_size_per_group)) {
    stop("Error: `sample_size_per_group` must be numeric with no NA.")
  }

  if (length(sample_size_per_group) == 1L) {
    if (sample_size_per_group <= 0) {
      stop("Error: `sample_size_per_group` must be > 0.")
    }

    n_by_arm <- rep(as.integer(sample_size_per_group), K)

  } else if (length(sample_size_per_group) == K) {
    if (any(sample_size_per_group <= 0)) {
      stop("Error: all `sample_size_per_group` entries must be > 0.")
    }

    n_by_arm <- as.integer(sample_size_per_group)

  } else {
    stop(
      "Error: `sample_size_per_group` must have length 1 or length K ",
      "(K = ", K, ")."
    )
  }

  # ---- per-endpoint checks -------------------------------------------------
  for (j in seq_along(endpoint_details)) {
    spec <- endpoint_details[[j]]
    typ  <- endpoint_types[j]

    if (typ == "binary" &&
        !is.null(spec$trt_effect) &&
        !is.null(spec$trt_prob)) {
      stop(
        "Error: Binary endpoint ", j, ": specify only one of ",
        "`trt_effect` or `trt_prob`."
      )
    }

    if (typ == "count" &&
        !is.null(spec$trt_effect) &&
        !is.null(spec$trt_count)) {
      stop(
        "Error: Count endpoint ", j, ": specify only one of ",
        "`trt_effect` or `trt_count`."
      )
    }

    if (typ == "continuous") {
      if (is.null(spec$baseline_mean) || is.null(spec$sd)) {
        stop("Error: Continuous endpoint ", j, " requires `baseline_mean` and `sd`.")
      }

      if (!is.numeric(spec$baseline_mean) || length(spec$baseline_mean) != 1L ||
          is.na(spec$baseline_mean)) {
        stop(
          "Error: Continuous endpoint ", j,
          " `baseline_mean` must be a single numeric value."
        )
      }

      if (!is.numeric(spec$sd) || anyNA(spec$sd)) {
        stop("Error: Continuous endpoint ", j, " `sd` must be numeric with no NA.")
      }

      if (!(length(spec$sd) %in% c(1L, K))) {
        stop(
          "Error: Continuous endpoint ", j, " `sd` must have length 1 or K ",
          "(K = ", K, ")."
        )
      }

      if (any(spec$sd <= 0)) {
        stop("Error: Continuous endpoint ", j, " `sd` must be > 0.")
      }

    } else if (typ == "binary") {
      if (is.null(spec$baseline_prob)) {
        stop("Error: Binary endpoint ", j, " requires `baseline_prob`.")
      }

      if (!is.numeric(spec$baseline_prob) || length(spec$baseline_prob) != 1L ||
          is.na(spec$baseline_prob)) {
        stop(
          "Error: Binary endpoint ", j,
          " `baseline_prob` must be a single numeric value."
        )
      }

      if (!(spec$baseline_prob > 0 && spec$baseline_prob < 1)) {
        stop("Error: Binary endpoint ", j, " `baseline_prob` must be in (0, 1).")
      }

      if (!is.null(spec$trt_prob)) {
        if (!is.numeric(spec$trt_prob) || anyNA(spec$trt_prob)) {
          stop("Error: Binary endpoint ", j, " `trt_prob` must be numeric with no NA.")
        }

        if (!(length(spec$trt_prob) %in% c(1L, K - 1L))) {
          stop(
            "Error: Binary endpoint ", j,
            " `trt_prob` must have length 1 or K-1 ",
            "(K = ", K, ")."
          )
        }

        if (any(spec$trt_prob <= 0 | spec$trt_prob >= 1)) {
          stop("Error: Binary endpoint ", j, " `trt_prob` must be in (0, 1).")
        }
      }

    } else if (typ == "count") {
      if (is.null(spec$baseline_mean) || is.null(spec$size)) {
        stop("Error: Count endpoint ", j, " requires `baseline_mean` and `size`.")
      }

      if (!is.numeric(spec$baseline_mean) || length(spec$baseline_mean) != 1L ||
          is.na(spec$baseline_mean) || spec$baseline_mean <= 0) {
        stop(
          "Error: Count endpoint ", j,
          " `baseline_mean` must be a single positive numeric value."
        )
      }

      if (!is.numeric(spec$size) || length(spec$size) != 1L ||
          is.na(spec$size) || spec$size <= 0) {
        stop(
          "Error: Count endpoint ", j,
          " `size` must be a single positive numeric value."
        )
      }

      if (!is.null(spec$p_zero)) {
        if (!is.numeric(spec$p_zero) || length(spec$p_zero) != 1L ||
            is.na(spec$p_zero)) {
          stop(
            "Error: Count endpoint ", j,
            " `p_zero` must be a single numeric value."
          )
        }

        if (spec$p_zero < 0 || spec$p_zero > 1) {
          stop("Error: Count endpoint ", j, " `p_zero` must be in [0, 1].")
        }
      }

      if (!is.null(spec$trt_count)) {
        if (!is.numeric(spec$trt_count) || anyNA(spec$trt_count)) {
          stop("Error: Count endpoint ", j, " `trt_count` must be numeric with no NA.")
        }

        if (!(length(spec$trt_count) %in% c(1L, K - 1L))) {
          stop(
            "Error: Count endpoint ", j,
            " `trt_count` must have length 1 or K-1 ",
            "(K = ", K, ")."
          )
        }

        if (any(spec$trt_count <= 0)) {
          stop("Error: Count endpoint ", j, " `trt_count` must be > 0.")
        }
      }

    } else if (typ == "time-to-event") {
      if (is.null(spec$baseline_rate)) {
        stop("Error: TTE endpoint ", j, " requires `baseline_rate`.")
      }

      if (!is.numeric(spec$baseline_rate) || length(spec$baseline_rate) != 1L ||
          is.na(spec$baseline_rate) || spec$baseline_rate <= 0) {
        stop(
          "Error: TTE endpoint ", j,
          " `baseline_rate` must be a single positive numeric value."
        )
      }

      if (!is.null(spec$censoring_rate)) {
        if (!is.numeric(spec$censoring_rate) ||
            length(spec$censoring_rate) != 1L ||
            is.na(spec$censoring_rate) ||
            spec$censoring_rate < 0) {
          stop(
            "Error: TTE endpoint ", j,
            " `censoring_rate` must be a single numeric value >= 0."
          )
        }
      }

      if (!is.null(spec$fatal_event)) {
        if (!is.logical(spec$fatal_event) || length(spec$fatal_event) != 1L ||
            is.na(spec$fatal_event)) {
          stop("Error: TTE endpoint ", j, " `fatal_event` must be TRUE or FALSE.")
        }
      }
    }
  }

  # ---- fatal event structure checks ----------------------------------------
  tte_idx <- which(endpoint_types == "time-to-event")

  if (length(tte_idx) > 0L) {
    fatal_flags <- vapply(
      tte_idx,
      function(j) isTRUE(endpoint_details[[j]]$fatal_event %||% FALSE),
      logical(1)
    )

    if (sum(fatal_flags) > 2L) {
      stop("Error: Only two terminal/fatal events are allowed.")
    }

    if (sum(fatal_flags) > 0L) {
      fatal_pos <- which(fatal_flags)

      if (any(fatal_pos > 2L)) {
        stop(
          "Error: Terminal/fatal events must be the first or first and second ",
          "TTE endpoints."
        )
      }
    }

    if (any(diff(as.integer(fatal_flags)) == 1L)) {
      stop(
        "Error: Fatal TTE endpoints must be listed before any non-fatal TTE ",
        "endpoints among TTE endpoints."
      )
    }
  }

  # ---- non_fatal_censors_fatal ---------------------------------------------
  if (!is.logical(non_fatal_censors_fatal) ||
      length(non_fatal_censors_fatal) != 1L ||
      is.na(non_fatal_censors_fatal)) {
    stop("Error: `non_fatal_censors_fatal` must be TRUE or FALSE.")
  }

  # ---- enrollment details --------------------------------------------------
  if (!is.list(enrollment_details)) {
    stop("Error: `enrollment_details` must be a list.")
  }

  dist <- enrollment_details$enrollment_distribution %||% "none"

  valid_dists <- c("none", "exponential", "piecewise")

  if (!dist %in% valid_dists) {
    stop(
      "Error: `enrollment_distribution` must be one of: ",
      paste(valid_dists, collapse = ", "),
      "."
    )
  }

  if (dist == "exponential") {
    rate <- enrollment_details$enrollment_exponential_rate %||% NULL

    if (is.null(rate)) {
      stop(
        "Error: `enrollment_exponential_rate` must be specified when ",
        "`enrollment_distribution = 'exponential'`."
      )
    }

    if (!is.numeric(rate) || length(rate) != 1L ||
        is.na(rate) || rate <= 0) {
      stop("Error: `enrollment_exponential_rate` must be a single positive number.")
    }
  }

  if (dist == "piecewise") {
    cuts  <- enrollment_details$piecewise_enrollment_cutpoints %||% NULL
    rates <- enrollment_details$piecewise_enrollment_rates %||% NULL

    if (is.null(cuts) || is.null(rates)) {
      stop(
        "Error: `piecewise_enrollment_cutpoints` and ",
        "`piecewise_enrollment_rates` must both be specified when ",
        "`enrollment_distribution = 'piecewise'`."
      )
    }

    if (!is.numeric(cuts) || length(cuts) < 2L || anyNA(cuts)) {
      stop(
        "Error: `piecewise_enrollment_cutpoints` must be numeric with at ",
        "least two values and no NA."
      )
    }

    if (cuts[1] != 0) {
      stop("Error: `piecewise_enrollment_cutpoints` should start at 0.")
    }

    if (any(diff(cuts) <= 0)) {
      stop("Error: `piecewise_enrollment_cutpoints` must be strictly increasing.")
    }

    if (!is.numeric(rates) || anyNA(rates)) {
      stop("Error: `piecewise_enrollment_rates` must be numeric with no NA.")
    }

    if (length(rates) != length(cuts) - 1L) {
      stop(
        "Error: `piecewise_enrollment_rates` must have length ",
        "`length(piecewise_enrollment_cutpoints) - 1`."
      )
    }

    if (any(rates < 0)) {
      stop("Error: `piecewise_enrollment_rates` must be non-negative.")
    }

    if (all(rates == 0)) {
      stop("Error: At least one `piecewise_enrollment_rates` value must be > 0.")
    }

    if (tail(rates, 1L) <= 0) {
      warning(
        "The final piecewise enrollment rate is 0. If the finite piecewise ",
        "accrual window does not generate enough subjects, simulation will fail."
      )
    }
  }

  # ---- follow-up details ---------------------------------------------------
  if (!is.list(followup_details)) {
    stop("Error: `followup_details` must be a list.")
  }

  min_followup <- followup_details$min_followup %||% NULL
  max_followup <- followup_details$max_followup %||% NULL

  if (!is.null(min_followup)) {
    if (!is.numeric(min_followup) || length(min_followup) != 1L ||
        is.na(min_followup) || min_followup < 0) {
      stop(
        "Error: `followup_details$min_followup` must be a non-negative ",
        "numeric scalar."
      )
    }
  }

  if (!is.null(max_followup)) {
    if (!is.numeric(max_followup) || length(max_followup) != 1L ||
        is.na(max_followup) || max_followup <= 0) {
      stop(
        "Error: `followup_details$max_followup` must be a positive numeric scalar."
      )
    }
  }

  if (!is.null(min_followup) && !is.null(max_followup) &&
      min_followup > max_followup) {
    warning(
      "`followup_details$min_followup` is greater than ",
      "`followup_details$max_followup`. Check that this is intentional."
    )
  }

  # ---- trial-end details ---------------------------------------------------
  if (!is.list(trial_end_details)) {
    stop("Error: `trial_end_details` must be a list.")
  }

  trial_end_type <- trial_end_details$type %||% "none"

  valid_trial_end_types <- c(
    "none",
    "fixed_calendar",
    "last_patient_min_followup",
    "event_driven"
  )

  if (!trial_end_type %in% valid_trial_end_types) {
    stop(
      "Error: `trial_end_details$type` must be one of: ",
      paste(valid_trial_end_types, collapse = ", "),
      "."
    )
  }

  if (trial_end_type == "fixed_calendar") {
    trial_end_time <- trial_end_details$trial_end_time %||% NULL

    if (is.null(trial_end_time) ||
        !is.numeric(trial_end_time) ||
        length(trial_end_time) != 1L ||
        is.na(trial_end_time) ||
        trial_end_time < 0) {
      stop(
        "Error: For `trial_end_details$type = 'fixed_calendar'`, ",
        "`trial_end_details$trial_end_time` must be a non-negative numeric scalar."
      )
    }
  }

  if (trial_end_type == "last_patient_min_followup") {
    if (is.null(min_followup)) {
      stop(
        "Error: `followup_details$min_followup` is required when ",
        "`trial_end_details$type = 'last_patient_min_followup'`."
      )
    }
  }

  if (trial_end_type == "event_driven") {
    if (length(tte_idx) == 0L) {
      stop(
        "Error: `trial_end_details$type = 'event_driven'` requires at least ",
        "one time-to-event endpoint."
      )
    }

    event_endpoint <- trial_end_details$event_endpoint %||% NULL
    target_events  <- trial_end_details$target_events %||% NULL

    # Validate event endpoint format/range.
    resolve_event_endpoint(event_endpoint = event_endpoint, tte_idx = tte_idx)

    if (is.null(target_events) ||
        !is.numeric(target_events) ||
        length(target_events) != 1L ||
        is.na(target_events) ||
        target_events < 1) {
      stop(
        "Error: `trial_end_details$target_events` must be a positive numeric ",
        "scalar when `trial_end_details$type = 'event_driven'`."
      )
    }

    if (target_events != as.integer(target_events)) {
      warning("`trial_end_details$target_events` will be coerced to integer.")
    }

    require_min_followup <- trial_end_details$require_min_followup %||% FALSE

    if (!is.logical(require_min_followup) ||
        length(require_min_followup) != 1L ||
        is.na(require_min_followup)) {
      stop("Error: `trial_end_details$require_min_followup` must be TRUE or FALSE.")
    }

    if (isTRUE(require_min_followup) && is.null(min_followup)) {
      stop(
        "Error: `trial_end_details$require_min_followup = TRUE` requires ",
        "`followup_details$min_followup`."
      )
    }

    max_trial_duration <- trial_end_details$max_trial_duration %||% NULL

    if (!is.null(max_trial_duration)) {
      if (!is.numeric(max_trial_duration) ||
          length(max_trial_duration) != 1L ||
          is.na(max_trial_duration) ||
          max_trial_duration <= 0) {
        stop(
          "Error: `trial_end_details$max_trial_duration` must be a positive ",
          "numeric scalar."
        )
      }
    }

    target_not_reached <- trial_end_details$target_not_reached %||% "error"

    valid_target_not_reached <- c(
      "error",
      "max_trial_duration",
      "last_patient_min_followup"
    )

    if (!target_not_reached %in% valid_target_not_reached) {
      stop(
        "Error: `trial_end_details$target_not_reached` must be one of: ",
        paste(valid_target_not_reached, collapse = ", "),
        "."
      )
    }

    if (target_not_reached == "max_trial_duration" &&
        is.null(max_trial_duration)) {
      stop(
        "Error: `target_not_reached = 'max_trial_duration'` requires ",
        "`trial_end_details$max_trial_duration`."
      )
    }

    if (target_not_reached == "last_patient_min_followup" &&
        is.null(min_followup)) {
      stop(
        "Error: `target_not_reached = 'last_patient_min_followup'` requires ",
        "`followup_details$min_followup`."
      )
    }
  }

  # ---- target_correlation --------------------------------------------------
  if (!is.logical(target_correlation) || length(target_correlation) != 1L ||
      is.na(target_correlation)) {
    stop("Error: `target_correlation` must be TRUE or FALSE.")
  }

  # ---- calibration checks --------------------------------------------------
  if (isTRUE(target_correlation) && !single_endpoint_mode) {
    if (!is.list(calibration_control)) {
      stop("Error: `calibration_control` must be a list when target_correlation = TRUE.")
    }

    n_mc      <- calibration_control$n_mc      %||% calibration_control$n_obs      %||% 10000
    tol       <- calibration_control$tol       %||% 0.001
    maxit     <- calibration_control$maxit     %||% 100
    rho_cap   <- calibration_control$rho_cap   %||% calibration_control$rho_max    %||% 0.999
    ensure_pd <- calibration_control$ensure_pd %||% calibration_control$ensure_cor_mat %||% TRUE
    convT     <- calibration_control$conv_norm_type %||% "F"

    if (!is.numeric(n_mc) || length(n_mc) != 1L || is.na(n_mc) || n_mc <= 0) {
      stop("Error: `calibration_control$n_mc` or legacy `n_obs` must be > 0.")
    }

    if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
      stop("Error: `calibration_control$tol` must be > 0.")
    }

    if (!is.numeric(maxit) || length(maxit) != 1L || is.na(maxit) || maxit <= 0) {
      stop("Error: `calibration_control$maxit` must be > 0.")
    }

    if (!is.numeric(rho_cap) || length(rho_cap) != 1L ||
        is.na(rho_cap) || rho_cap <= 0 || rho_cap >= 1) {
      stop(
        "Error: `calibration_control$rho_cap` or legacy `rho_max` must be ",
        "in (0, 1)."
      )
    }

    if (!is.logical(ensure_pd) || length(ensure_pd) != 1L ||
        is.na(ensure_pd)) {
      stop(
        "Error: `calibration_control$ensure_pd` or legacy ",
        "`ensure_cor_mat` must be TRUE or FALSE."
      )
    }

    if (!is.character(convT) || length(convT) != 1L || is.na(convT)) {
      stop("Error: `calibration_control$conv_norm_type` must be a single character value.")
    }
  }

  invisible(list(
    n_arms = K,
    n_by_arm = n_by_arm,
    endpoint_types = endpoint_types,
    control_only = control_only,
    single_endpoint_mode = single_endpoint_mode
  ))
}
