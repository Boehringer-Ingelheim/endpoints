##############################################################################
# Blinded estimate of placebo 1-year event rate from pooled KM event rate
##############################################################################
estimate_blinded_pbo_1yr_event_rate_and_SSR <- function(
    dat,
    time_col = "TTE_1",
    status_col = "Status_1",

    # Planning assumptions used to decode the blinded pooled rate
    assumed_pbo_2yr_event_rate = 0.65,
    assumed_hr = 0.55,
    allocation_ratio = c(3, 2),

    # Decision rule
    target_poe = 0.85,
    ssr_pbo_1yr_cutoff = 1 - sqrt(1 - 0.65),

    timepoint_months = 12,
    tolerance = 1e-8,
    decoding_method = c(
      "km",
      "exponential",
      "km_extrapolate"
    )
) {

  decoding_method <- match.arg(decoding_method)

  # Planned placebo 1-year event rate under exponential survival.
  # If p_2yr = 0.65, then S(2) = 1 - 0.65 = 0.35,
  # and S(1) = sqrt(S(2)).
  assumed_pbo_1yr_event_rate <- 1 - sqrt(1 - assumed_pbo_2yr_event_rate)

  # Extract once for readability
  times    <- dat[[time_col]]
  statuses <- dat[[status_col]]

  n_events_total    <- sum(statuses, na.rm = TRUE)
  total_person_time <- sum(times, na.rm = TRUE)

  # ------------------------------------------------------------
  # 1. Estimate pooled 1-year event rate by chosen method
  # ------------------------------------------------------------
  if (decoding_method == "km") {

    # Kaplan-Meier estimate of S(timepoint_months), with extend = TRUE
    # Pros: nonparametric, robust to non-exponential hazards.
    # Cons: when BSSR fires before most patients reach `timepoint_months`
    #       of follow-up, extend = TRUE flatlines S(timepoint_months) and
    #       biases the 1-year event rate downward.
    if (n_events_total == 0L) {
      pooled_1yr_event_rate <- 0
    } else {
      fit <- survival::survfit(
        survival::Surv(time = times, event = statuses) ~ 1
      )
      surv_at_timepoint <- summary(
        fit,
        times  = timepoint_months,
        extend = TRUE
      )$surv
      pooled_1yr_event_rate <- 1 - surv_at_timepoint
    }

  } else if (decoding_method == "km_extrapolate") { # this is new

    # Anchored KM (9.75 months) + constant-hazard extrapolation to timepoint.
    #
    # At an early BSSR the KM curve is well supported out to ~9.75 months but is
    # flatlined beyond the last observed event if read directly at 12 months.
    # Instead we read KM at 9.75 months, convert to an average monthly hazard over
    # [0, 9.75], and project that hazard forward to the timepoint:
    #
    #   S_hat(6)  from KM
    #   lambda_bar = -log(S_hat(6)) / 6            (avg monthly hazard)
    #   S_hat(12)  = exp(-lambda_bar * timepoint_months)
    #            = S_hat(6)^(timepoint_months / 6)
    #
    # Pros: uses the well-estimated early portion of the curve and projects
    #       smoothly, avoiding the downward bias from flatlining a sparse tail.
    # Cons: assumes (locally) constant hazard between 6 months and the timepoint.
    anchor_months <- 9.75

    if (n_events_total == 0L) {
      pooled_1yr_event_rate <- 0
    } else {
      fit <- survival::survfit(
        survival::Surv(time = times, event = statuses) ~ 1
      )
      surv_at_anchor <- summary(
        fit,
        times  = anchor_months,
        extend = TRUE
      )$surv

      # Guard against degenerate S(ancho_month) values
      if (is.na(surv_at_anchor) || surv_at_anchor >= 1) {
        # No events observed by anchor time -> no information to project a rate
        pooled_1yr_event_rate <- 0
      } else if (surv_at_anchor <= 0) {
        pooled_1yr_event_rate <- 1
      } else {
        # Project survival from 6 months to the timepoint under constant hazard
        surv_at_timepoint <- surv_at_anchor^(timepoint_months / anchor_months)
        pooled_1yr_event_rate <- 1 - surv_at_timepoint
      }
    }} else if (decoding_method == "exponential") {

      # Pooled exponential MLE:
      #   lambda_hat = total events / total person-time
      #   p_1yr_hat  = 1 - exp(-lambda_hat * timepoint_months)
      if (n_events_total == 0L || total_person_time <= 0) {
        pooled_1yr_event_rate <- 0
      } else {
        lambda_hat            <- n_events_total / total_person_time
        pooled_1yr_event_rate <- 1 - exp(-lambda_hat * timepoint_months)
      }
    }

  # ------------------------------------------------------------
  # 2. Decode pooled rate to implied placebo rate
  #
  # Under blinding:
  #
  #   p_pool = w_trt * p_trt + w_pbo * p_pbo
  #
  # Exponential model:
  #
  #   p_pbo = 1 - exp(-lambda_pbo)
  #   p_trt = 1 - exp(-lambda_trt)
  #   lambda_trt = assumed_hr * lambda_pbo
  #
  # Solve numerically for p_pbo.
  # ------------------------------------------------------------

  w_trt <- allocation_ratio[1] / sum(allocation_ratio)
  w_pbo <- allocation_ratio[2] / sum(allocation_ratio)

  encode_pooled_rate <- function(pbo_1yr_event_rate) {

    lambda_pbo <- -log(1 - pbo_1yr_event_rate)
    lambda_trt <- assumed_hr * lambda_pbo

    p_pbo <- 1 - exp(-lambda_pbo)
    p_trt <- 1 - exp(-lambda_trt)

    w_trt * p_trt + w_pbo * p_pbo
  }

  # these lines catch bad cases, but shouldn't be needed
  if (is.na(pooled_1yr_event_rate)) {
    pbo_1yr_event_rate_hat <- NA_real_
  } else if (pooled_1yr_event_rate <= 0) {
    pbo_1yr_event_rate_hat <- 0
  } else if (pooled_1yr_event_rate >= 1) {
    pbo_1yr_event_rate_hat <- 1
  } else {
    lower <- 0
    upper <- 1
    # small algorithm to numerically find the implied PBO rate from the data
    while ((upper - lower) > tolerance) {
      mid        <- (lower + upper) / 2
      pooled_mid <- encode_pooled_rate(mid)
      if (pooled_mid > pooled_1yr_event_rate) {
        upper <- mid
      } else {
        lower <- mid
      }
    }
    pbo_1yr_event_rate_hat <- (lower + upper) / 2
  }

  trigger_ssr <-
    pbo_1yr_event_rate_hat < ssr_pbo_1yr_cutoff
  # from the protocol, this read 'If the projected probability of achieving efficacy
  # by the end of the trial is below 85%, the sample size will be increased.' However,
  # this was essentially assessed in the first set of simulations
  # drop_fpf_recruitment.R, enroll_sim_endpoints_noMonthRounding.R, drop_sim_endpoints.R
  # in those sims, efficacy fell below 85% if the encoded PBO event rate was low
  # here, we assume a HR=.55 and proportional hazards, same prior etc, so the sims match

  tibble(
    n_blinded = nrow(dat),
    n_events_total = n_events_total,
    total_person_time = total_person_time,
    timepoint_months = timepoint_months,

    decoding_method = decoding_method,
    pooled_1yr_event_rate = pooled_1yr_event_rate,

    assumed_hr = assumed_hr,
    assumed_pbo_2yr_event_rate =
      assumed_pbo_2yr_event_rate,
    assumed_pbo_1yr_event_rate =
      1 - sqrt(1 - assumed_pbo_2yr_event_rate),

    target_poe = target_poe,
    ssr_pbo_1yr_cutoff = ssr_pbo_1yr_cutoff,
    ssr_pbo_2yr_cutoff =
      1 - (1 - ssr_pbo_1yr_cutoff)^2,

    w_trt = w_trt,
    w_pbo = w_pbo,

    pbo_1yr_event_rate_hat =
      pbo_1yr_event_rate_hat,

    pbo_2yr_event_rate_hat =
      1 - (1 - pbo_1yr_event_rate_hat)^2,

    trigger_ssr = trigger_ssr
  )
}
