# purpose & overview
#
# 1) "marginals and correlation structure look correct (large n, multi-arm)"
#    Generates a mixed-endpoint, multi-arm trial in a large sample and checks that:
#      - observed within-arm Pearson correlation matrices are close to the
#        user-supplied target correlation matrix,
#      - each endpoint’s marginal distribution matches its inputs.
#
# 2) "TTE censoring_rate supports scalar broadcast and arm-specific vectors"
#    Checks that a scalar censoring rate is used for all arms and a vector
#    censoring rate is applied arm-by-arm.
#
# 3) "semi-competing risks: non-fatal cannot occur after fatal (multi-arm)"
#    Validates the semi-competing risks rule: if a non-fatal endpoint is observed,
#    it must occur on or before the fatal time.
#
# 4) "non_fatal_censors_fatal=TRUE: censoring of non-fatal censors other TTE endpoints (multi-arm)"
#    Validates the optional rule that censoring of a non-fatal endpoint can censor
#    other TTE endpoints.
#
# 5) "argument checks: required fields and structural rules enforced"
#    Ensures validation is working for required fields, structural TTE rules,
#    arm-length rules, and mutually exclusive treatment inputs.
#
# 6) "argument checks: TTE censoring_rate length and values enforced"
#
# 7) "trial calendar: fixed calendar uses enrollTime-adjusted available follow-up"
#    Validates stochastic enrollment + fixed calendar trial end:
#      availableFollowup = pmin(max_followup, pmax(0, trial_end_time - enrollTime)).
#
# 8) "trial calendar: max_followup cap binds when trial end is late"
#    Validates that max_followup caps individual follow-up.
#
# 9) "enrollment and trial-calendar argument dependencies are enforced"
#    Ensures new enrollment/follow-up/trial-end settings are internally consistent.
#
# 10) "event-driven trial end: target events reached and analysis population identified"
#    Validates event-driven trial end and checks use of availableFollowup > 0.
#
# 11) "event-driven trial end: require_min_followup delays final analysis"
#    Validates event-driven plus last-patient minimum follow-up.
#
# 12) "event-driven trial end: target-not-reached fallback works"
#     Validates fallback to last-patient-min-follow-up.
#
# 13) "arm_mode='control': omit trt column and behave like K=1"
#
# 14) "single-endpoint mode: correlation_matrix=NULL generates one endpoint"
# -------------------------------------------------------------------------


# --- helpers ---------------------------------------------------------------

skip_if_missing_pkgs <- function() {
  testthat::skip_if_not_installed("MASS")
  testthat::skip_if_not_installed("survival")
}

make_cor4 <- function() {
  corr_make(
    num_endpoints = 4,
    values = rbind(
      c(1, 2, 0.20),
      c(1, 3, 0.10),
      c(1, 4, 0.15),
      c(2, 3, 0.25),
      c(2, 4, 0.05),
      c(3, 4, 0.30)
    )
  )
}

offdiag_max_abs <- function(A, B) {
  idx <- which(row(A) != col(A), arr.ind = TRUE)
  max(abs(A[idx] - B[idx]))
}

cal_ctl_fast <- list(
  n_mc = 5000,
  tol = 0.005,
  maxit = 75,
  rho_cap = 0.999,
  ensure_pd = TRUE,
  conv_norm_type = "O"
)

get_cor_arm <- function(s, k) {
  nm <- paste0("arm_", k)
  if (!is.null(s$estimated_correlation_by_arm) &&
      !is.null(s$estimated_correlation_by_arm[[nm]])) {
    return(s$estimated_correlation_by_arm[[nm]])
  }
  NULL
}

ep_tte_calendar <- list(
  endpoint_type  = "tte",
  baseline_rate  = 1 / 52,
  trt_effect     = log(0.8),
  censoring_rate = 0,
  fatal_event    = TRUE
)

ep_tte_fast <- list(
  endpoint_type  = "tte",
  baseline_rate  = 1 / 12,
  trt_effect     = log(0.8),
  censoring_rate = 0,
  fatal_event    = TRUE
)


# --- large-sample marginals + correlation ----------------------------------

testthat::test_that("marginals and correlation structure look correct (large n, multi-arm)", {
  skip_if_missing_pkgs()

  set.seed(1)
  cor_mat <- make_cor4()

  n_by_arm <- c(5000, 6500, 5500)
  sd_vec <- c(3.0, 4.0, 2.5)

  obj <- makeData(
    correlation_matrix = cor_mat,
    sample_size_per_group = n_by_arm,
    SEED = 123,
    arm_mode = "full",
    endpoint_details = list(
      list(
        endpoint_type = "normal",
        baseline_mean = 10,
        sd = sd_vec,
        trt_effect = c(-2, -1)
      ),
      list(
        endpoint_type = "binary",
        baseline_prob = 0.30,
        trt_prob = c(0.25, 0.20)
      ),
      list(
        endpoint_type = "count",
        baseline_mean = 8,
        size = 5,
        p_zero = 0,
        trt_count = c(10, 9)
      ),
      list(
        endpoint_type = "tte",
        baseline_rate = 1 / 24,
        censoring_rate = 0,
        fatal_event = FALSE,
        trt_effect = c(log(0.80), log(0.90))
      )
    ),
    enrollment_details = list(
      enrollment_distribution = "none"
    ),
    followup_details = list(),
    trial_end_details = list(
      type = "none"
    ),
    non_fatal_censors_fatal = FALSE,
    target_correlation = TRUE,
    calibration_control = cal_ctl_fast
  )

  s <- summary(obj)
  d <- obj$data

  testthat::expect_true("trt" %in% names(d))

  tab <- table(d$trt)
  testthat::expect_equal(as.integer(tab["0"]), n_by_arm[1])
  testthat::expect_equal(as.integer(tab["1"]), n_by_arm[2])
  testthat::expect_equal(as.integer(tab["2"]), n_by_arm[3])

  for (k in 0:2) {
    cor_k <- get_cor_arm(s, k)
    testthat::expect_true(!is.null(cor_k))
    testthat::expect_true(offdiag_max_abs(cor_k, cor_mat) < 0.07)
  }

  testthat::expect_true(abs(mean(d$Cont_1[d$trt == 0]) - 10) < 0.25)
  testthat::expect_true(abs(mean(d$Cont_1[d$trt == 1]) -  8) < 0.30)
  testthat::expect_true(abs(mean(d$Cont_1[d$trt == 2]) -  9) < 0.30)

  testthat::expect_true(abs(sd(d$Cont_1[d$trt == 0]) - sd_vec[1]) < 0.30)
  testthat::expect_true(abs(sd(d$Cont_1[d$trt == 1]) - sd_vec[2]) < 0.35)
  testthat::expect_true(abs(sd(d$Cont_1[d$trt == 2]) - sd_vec[3]) < 0.30)

  p0 <- mean(d$Bin_1[d$trt == 0] == 1)
  p1 <- mean(d$Bin_1[d$trt == 1] == 1)
  p2 <- mean(d$Bin_1[d$trt == 2] == 1)

  testthat::expect_true(abs(p0 - 0.30) < 0.02)
  testthat::expect_true(abs(p1 - 0.25) < 0.03)
  testthat::expect_true(abs(p2 - 0.20) < 0.03)

  m0 <- mean(d$Int_1[d$trt == 0])
  m1 <- mean(d$Int_1[d$trt == 1])
  m2 <- mean(d$Int_1[d$trt == 2])

  testthat::expect_true(abs(m0 -  8) < 0.40)
  testthat::expect_true(abs(m1 - 10) < 0.55)
  testthat::expect_true(abs(m2 -  9) < 0.55)

  t0 <- mean(d$TTE_1[d$trt == 0])
  t1 <- mean(d$TTE_1[d$trt == 1])
  t2 <- mean(d$TTE_1[d$trt == 2])

  testthat::expect_true(t1 > t0)
  testthat::expect_true(t2 > t0)
})


# --- arm-specific TTE censoring ---------------------------------------------

testthat::test_that("TTE censoring_rate supports scalar broadcast and arm-specific vectors", {
  baseline_rate <- 0.1
  n_by_arm <- c(5000, 5000, 5000)

  obj_vec <- makeData(
    correlation_matrix = NULL,
    sample_size_per_group = n_by_arm,
    SEED = 321,
    endpoint_details = list(
      list(
        endpoint_type = "tte",
        baseline_rate = baseline_rate,
        censoring_rate = c(0, 0.15, 0.6),
        fatal_event = FALSE
      )
    ),
    enrollment_details = list(
      enrollment_distribution = "none"
    ),
    followup_details = list(),
    trial_end_details = list(
      type = "none"
    ),
    non_fatal_censors_fatal = FALSE,
    target_correlation = FALSE
  )

  d_vec <- obj_vec$data
  obs_event_vec <- as.numeric(tapply(d_vec$Status_1, d_vec$trt, mean))
  target_event_vec <- baseline_rate / (baseline_rate + c(0, 0.15, 0.6))

  testthat::expect_equal(obj_vec$meta$n_arms, 3L)
  testthat::expect_equal(as.integer(table(d_vec$trt)), n_by_arm)
  testthat::expect_equal(obs_event_vec[1], 1)
  testthat::expect_true(all(diff(obs_event_vec) < 0))
  testthat::expect_true(max(abs(obs_event_vec[-1] - target_event_vec[-1])) < 0.04)

  scalar_censoring_rate <- 0.2
  obj_scalar <- makeData(
    correlation_matrix = NULL,
    sample_size_per_group = n_by_arm,
    SEED = 654,
    arm_mode = "full",
    endpoint_details = list(
      list(
        endpoint_type = "tte",
        baseline_rate = baseline_rate,
        censoring_rate = scalar_censoring_rate,
        fatal_event = FALSE,
        trt_effect = c(0, 0)
      )
    ),
    enrollment_details = list(
      enrollment_distribution = "none"
    ),
    followup_details = list(),
    trial_end_details = list(
      type = "none"
    ),
    non_fatal_censors_fatal = FALSE,
    target_correlation = FALSE
  )

  d_scalar <- obj_scalar$data
  obs_event_scalar <- as.numeric(tapply(d_scalar$Status_1, d_scalar$trt, mean))
  target_event_scalar <- baseline_rate / (baseline_rate + scalar_censoring_rate)

  testthat::expect_true(max(abs(obs_event_scalar - target_event_scalar)) < 0.035)
  testthat::expect_true(diff(range(obs_event_scalar)) < 0.05)
})


# --- semi-competing risks ---------------------------------------------------

testthat::test_that("semi-competing risks: non-fatal cannot occur after fatal (multi-arm)", {
  skip_if_missing_pkgs()

  set.seed(2)
  cor_mat <- corr_make(num_endpoints = 2, values = rbind(c(1, 2, 0.25)))
  n_by_arm <- c(4000, 4000, 4000)

  obj <- makeData(
    correlation_matrix = cor_mat,
    sample_size_per_group = n_by_arm,
    SEED = 789,
    arm_mode = "full",
    endpoint_details = list(
      list(
        endpoint_type = "tte",
        baseline_rate = 1 / 20,
        censoring_rate = 1 / 80,
        fatal_event = TRUE,
        trt_effect = c(log(0.9), log(0.95))
      ),
      list(
        endpoint_type = "tte",
        baseline_rate = 1 / 10,
        censoring_rate = 1 / 80,
        fatal_event = FALSE,
        trt_effect = c(0, 0)
      )
    ),
    enrollment_details = list(
      enrollment_distribution = "none"
    ),
    followup_details = list(),
    trial_end_details = list(
      type = "none"
    ),
    non_fatal_censors_fatal = FALSE,
    target_correlation = FALSE
  )

  d <- obj$data

  idx_nf_obs <- which(d$Status_2 == 1)

  if (length(idx_nf_obs) > 0) {
    testthat::expect_true(all(d$TTE_2[idx_nf_obs] <= d$TTE_1[idx_nf_obs] + 1e-12))
  }

  idx_bad <- which(d$Status_2 == 1 & d$TTE_2 > d$TTE_1 + 1e-12)
  testthat::expect_length(idx_bad, 0)
})


testthat::test_that("non_fatal_censors_fatal=TRUE: censoring of non-fatal censors other TTE endpoints (multi-arm)", {
  skip_if_missing_pkgs()

  set.seed(3)
  cor_mat <- corr_make(num_endpoints = 2, values = rbind(c(1, 2, 0.1)))
  n_by_arm <- c(3000, 3000, 3000)

  obj <- makeData(
    correlation_matrix = cor_mat,
    sample_size_per_group = n_by_arm,
    SEED = 999,
    arm_mode = "full",
    endpoint_details = list(
      list(
        endpoint_type = "tte",
        baseline_rate = 1 / 30,
        censoring_rate = 1 / 12,
        fatal_event = TRUE,
        trt_effect = c(0, 0)
      ),
      list(
        endpoint_type = "tte",
        baseline_rate = 1 / 25,
        censoring_rate = 1 / 12,
        fatal_event = FALSE,
        trt_effect = c(0, 0)
      )
    ),
    enrollment_details = list(
      enrollment_distribution = "none"
    ),
    followup_details = list(),
    trial_end_details = list(
      type = "none"
    ),
    non_fatal_censors_fatal = TRUE,
    target_correlation = FALSE
  )

  d <- obj$data

  idx_cens_nf <- which(d$Status_2 == 0)

  if (length(idx_cens_nf) > 0) {
    testthat::expect_true(all(d$TTE_1[idx_cens_nf] <= d$TTE_2[idx_cens_nf] + 1e-12))
  }
})


# --- argument requirement logic ---------------------------------------------

testthat::test_that("argument checks: required fields and structural rules enforced", {

  testthat::expect_error(
    makeData(
      correlation_matrix = NULL,
      sample_size_per_group = 10,
      SEED = 1,
      endpoint_details = list(
        list(endpoint_type = "normal", baseline_mean = 0, sd = 1),
        list(endpoint_type = "binary", baseline_prob = 0.3)
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE,
      arm_mode = "control"
    ),
    "exactly one endpoint|single.*endpoint|correlation_matrix",
    ignore.case = TRUE
  )

  cor_mat <- diag(1)

  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = 10,
      SEED = 1,
      endpoint_details = list(
        list(endpoint_type = "normal", baseline_mean = 0, trt_effect = 0)
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE,
      arm_mode = "full"
    ),
    "requires.*sd|sd",
    ignore.case = TRUE
  )

  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = 10,
      SEED = 1,
      endpoint_details = list(
        list(endpoint_type = "binary", baseline_prob = 1.0, trt_effect = 0)
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE,
      arm_mode = "full"
    ),
    "baseline_prob",
    ignore.case = TRUE
  )

  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = 10,
      SEED = 1,
      endpoint_details = list(
        list(endpoint_type = "count", baseline_mean = 5, trt_effect = 0)
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE,
      arm_mode = "full"
    ),
    "requires.*size|size",
    ignore.case = TRUE
  )

  cor_mat2 <- diag(3)

  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat2,
      sample_size_per_group = c(10, 10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, fatal_event = FALSE, trt_effect = c(0, 0)),
        list(endpoint_type = "tte", baseline_rate = 0.1, fatal_event = FALSE, trt_effect = c(0, 0)),
        list(endpoint_type = "tte", baseline_rate = 0.1, fatal_event = TRUE,  trt_effect = c(0, 0))
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "Terminal|fatal|first",
    ignore.case = TRUE
  )

  testthat::expect_error(
    makeData(
      correlation_matrix = diag(1),
      sample_size_per_group = c(10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "binary", baseline_prob = 0.3, trt_effect = 0.1, trt_prob = 0.2)
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "trt_prob|trt_effect|only one|exclusive",
    ignore.case = TRUE
  )

  testthat::expect_error(
    makeData(
      correlation_matrix = diag(1),
      sample_size_per_group = c(10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "count", baseline_mean = 5, size = 2, trt_effect = 0.1, trt_count = 6)
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "trt_count|trt_effect|only one|exclusive",
    ignore.case = TRUE
  )
})


# --- arm-specific TTE censoring argument checks -----------------------------

testthat::test_that("argument checks: TTE censoring_rate length and values enforced", {
  testthat::expect_error(
    makeData(
      correlation_matrix = NULL,
      sample_size_per_group = c(10, 10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(
          endpoint_type = "tte",
          baseline_rate = 0.1,
          trt_effect = c(0, 0),
          censoring_rate = c(0.01, 0.02)
        )
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "censoring_rate.*length 1 or K",
    ignore.case = TRUE
  )

  testthat::expect_error(
    makeData(
      correlation_matrix = NULL,
      sample_size_per_group = c(10, 10, 10),
      SEED = 1,
      endpoint_details = list(
        list(
          endpoint_type = "tte",
          baseline_rate = 0.1,
          censoring_rate = c(0.01, NA, 0.03)
        )
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "censoring_rate.*numeric with no NA",
    ignore.case = TRUE
  )

  testthat::expect_error(
    makeData(
      correlation_matrix = NULL,
      sample_size_per_group = c(10, 10, 10),
      SEED = 1,
      endpoint_details = list(
        list(
          endpoint_type = "tte",
          baseline_rate = 0.1,
          censoring_rate = c(0.01, -0.02, 0.03)
        )
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "censoring_rate.*>= 0",
    ignore.case = TRUE
  )

  testthat::expect_error(
    makeData(
      correlation_matrix = NULL,
      sample_size_per_group = 10,
      SEED = 1,
      arm_mode = "control",
      endpoint_details = list(
        list(
          endpoint_type = "tte",
          baseline_rate = 0.1,
          censoring_rate = c(0.01, 0.02)
        )
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "control-only.*censoring_rate",
    ignore.case = TRUE
  )
})


testthat::test_that("argument checks: sd length rules + sample_size length rules enforced", {

  cor_mat <- diag(1)

  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 12, 11),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(
          endpoint_type = "normal",
          baseline_mean = 0,
          sd = c(1, 1.2, 0.8),
          trt_effect = c(0.1, 0.2)
        )
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    NA
  )

  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 12, 11),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(
          endpoint_type = "normal",
          baseline_mean = 0,
          sd = c(1, 2),
          trt_effect = c(0.1, 0.2)
        )
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "sd.*length 1 or K|sd.*length",
    ignore.case = TRUE
  )

  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 12),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(
          endpoint_type = "normal",
          baseline_mean = 0,
          sd = 1,
          trt_effect = c(0.1, 0.2)
        )
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "sample_size_per_group|length",
    ignore.case = TRUE
  )
})


# --- enrollment / trial-calendar behavior -----------------------------------

testthat::test_that("trial calendar: fixed calendar uses enrollTime-adjusted available follow-up", {
  skip_if_missing_pkgs()

  cor_mat <- diag(1)
  trial_end_time <- 80
  max_followup <- 104

  obj <- makeData(
    correlation_matrix = cor_mat,
    sample_size_per_group = c(500, 500, 500),
    SEED = 88,
    arm_mode = "full",
    endpoint_details = list(
      list(
        endpoint_type = "tte",
        baseline_rate = 1 / 6,
        censoring_rate = 0,
        fatal_event = FALSE,
        trt_effect = c(0, 0)
      )
    ),
    enrollment_details = list(
      enrollment_distribution = "exponential",
      enrollment_exponential_rate = 20
    ),
    followup_details = list(
      max_followup = max_followup
    ),
    trial_end_details = list(
      type = "fixed_calendar",
      trial_end_time = trial_end_time
    ),
    non_fatal_censors_fatal = FALSE,
    target_correlation = FALSE
  )

  d <- obj$data
  tc <- obj$meta$trial_calendar

  testthat::expect_true("enrollTime" %in% names(d))
  testthat::expect_true("availableFollowup" %in% names(d))

  testthat::expect_equal(tc$trial_end_reason, "fixed_calendar")
  testthat::expect_equal(tc$trial_end_time, trial_end_time)

  testthat::expect_true(all(d$enrollTime >= 0))

  expected_fu <- pmin(max_followup, pmax(0, trial_end_time - d$enrollTime))

  testthat::expect_equal(d$availableFollowup, expected_fu, tolerance = 1e-10)
  testthat::expect_true(all(d$TTE_1 <= d$availableFollowup + 1e-12))

  idx0 <- which(d$availableFollowup == 0)

  if (length(idx0) > 0) {
    testthat::expect_true(all(abs(d$TTE_1[idx0]) < 1e-12))
    testthat::expect_true(all(d$Status_1[idx0] == 0))
  }

  idx_admin <- which(d$Status_1 == 0 & abs(d$TTE_1 - d$availableFollowup) < 1e-10)

  if (length(idx_admin) > 0) {
    testthat::expect_true(all(d$TTE_1[idx_admin] == d$availableFollowup[idx_admin]))
  }
})


testthat::test_that("trial calendar: max_followup cap binds when trial end is late", {
  skip_if_missing_pkgs()

  cor_mat <- diag(1)
  trial_end_time <- 200
  max_followup <- 104

  obj <- makeData(
    correlation_matrix = cor_mat,
    sample_size_per_group = c(400, 400, 400),
    SEED = 89,
    arm_mode = "full",
    endpoint_details = list(
      list(
        endpoint_type = "tte",
        baseline_rate = 1 / 40,
        censoring_rate = 0,
        fatal_event = FALSE,
        trt_effect = c(0, 0)
      )
    ),
    enrollment_details = list(
      enrollment_distribution = "exponential",
      enrollment_exponential_rate = 20
    ),
    followup_details = list(
      max_followup = max_followup
    ),
    trial_end_details = list(
      type = "fixed_calendar",
      trial_end_time = trial_end_time
    ),
    non_fatal_censors_fatal = FALSE,
    target_correlation = FALSE
  )

  d <- obj$data

  expected_fu <- pmin(max_followup, pmax(0, trial_end_time - d$enrollTime))

  testthat::expect_equal(d$availableFollowup, expected_fu, tolerance = 1e-10)
  testthat::expect_true(any(abs(d$availableFollowup - max_followup) < 1e-10))
  testthat::expect_true(all(d$availableFollowup <= max_followup + 1e-12))
  testthat::expect_true(all(d$TTE_1 <= d$availableFollowup + 1e-12))
})


testthat::test_that("trial calendar: last_patient_min_followup rule is applied", {
  skip_if_missing_pkgs()

  cor_mat <- diag(1)
  min_followup <- 52
  max_followup <- 156

  obj <- makeData(
    correlation_matrix = cor_mat,
    sample_size_per_group = c(300, 300, 300),
    SEED = 90,
    arm_mode = "full",
    endpoint_details = list(
      list(
        endpoint_type = "tte",
        baseline_rate = 1 / 30,
        censoring_rate = 0,
        fatal_event = FALSE,
        trt_effect = c(0, 0)
      )
    ),
    enrollment_details = list(
      enrollment_distribution = "exponential",
      enrollment_exponential_rate = 20
    ),
    followup_details = list(
      min_followup = min_followup,
      max_followup = max_followup
    ),
    trial_end_details = list(
      type = "last_patient_min_followup"
    ),
    non_fatal_censors_fatal = FALSE,
    target_correlation = FALSE
  )

  d <- obj$data
  tc <- obj$meta$trial_calendar

  expected_trial_end <- max(d$enrollTime) + min_followup
  expected_fu <- pmin(max_followup, pmax(0, expected_trial_end - d$enrollTime))

  testthat::expect_equal(tc$trial_end_reason, "last_patient_min_followup")
  testthat::expect_equal(tc$trial_end_time, expected_trial_end, tolerance = 1e-10)
  testthat::expect_equal(d$availableFollowup, expected_fu, tolerance = 1e-10)
  testthat::expect_true(max(d$availableFollowup) <= max_followup + 1e-12)
  testthat::expect_true(min(d$availableFollowup) <= min_followup + 1e-12)
})


testthat::test_that("enrollment and trial-calendar argument dependencies are enforced", {
  cor_mat <- diag(1)

  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, trt_effect = c(0, 0))
      ),
      enrollment_details = list(
        administrative_censoring = 12
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "administrative_censoring|no longer supported|followup_details|trial_end_details",
    ignore.case = TRUE
  )

  # uniform no longer supported
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, trt_effect = c(0, 0))
      ),
      enrollment_details = list(
        enrollment_distribution = "uniform"
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "enrollment_distribution|none|exponential|piecewise",
    ignore.case = TRUE
  )

  # exponential enrollment requires positive rate
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, trt_effect = c(0, 0))
      ),
      enrollment_details = list(
        enrollment_distribution = "exponential",
        enrollment_exponential_rate = NULL
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "exponential_rate|enrollment_exponential_rate",
    ignore.case = TRUE
  )

  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, trt_effect = c(0, 0))
      ),
      enrollment_details = list(
        enrollment_distribution = "exponential",
        enrollment_exponential_rate = -1
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "positive|> 0",
    ignore.case = TRUE
  )

  # stochastic enrollment no longer requires admin censoring
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, trt_effect = c(0, 0))
      ),
      enrollment_details = list(
        enrollment_distribution = "exponential",
        enrollment_exponential_rate = 1
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    NA
  )

  # piecewise enrollment checks
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, trt_effect = c(0, 0))
      ),
      enrollment_details = list(
        enrollment_distribution = "piecewise",
        piecewise_enrollment_cutpoints = c(0, 5, 4),
        piecewise_enrollment_rates = c(0.1, 0.2)
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "increasing",
    ignore.case = TRUE
  )

  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, trt_effect = c(0, 0))
      ),
      enrollment_details = list(
        enrollment_distribution = "piecewise",
        piecewise_enrollment_cutpoints = c(0, 5, 10),
        piecewise_enrollment_rates = c(0.1)
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "length",
    ignore.case = TRUE
  )

  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, trt_effect = c(0, 0))
      ),
      enrollment_details = list(
        enrollment_distribution = "piecewise",
        piecewise_enrollment_cutpoints = c(0, 5, 10),
        piecewise_enrollment_rates = c(0, 0)
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "At least one|must be > 0",
    ignore.case = TRUE
  )

  testthat::expect_no_error({
    obj <- makeData(
      correlation_matrix = NULL,
      sample_size_per_group = c(40, 40),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "normal", baseline_mean = 0, sd = 1, trt_effect = 0)
      ),
      enrollment_details = list(
        enrollment_distribution = "piecewise",
        piecewise_enrollment_cutpoints = c(0, 4.5, 6),
        piecewise_enrollment_rates = c(24 / 4.5, 7, 8.9)
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    )

    d <- as.data.frame(obj)

    testthat::expect_equal(nrow(d), 80)
    testthat::expect_true("enrollTime" %in% names(d))
    testthat::expect_true(all(is.finite(d$enrollTime)))
    testthat::expect_equal(min(d$enrollTime), 0)
    testthat::expect_gt(max(d$enrollTime), 6)
  })

  # fixed calendar requires trial_end_time
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, trt_effect = 0)
      ),
      trial_end_details = list(
        type = "fixed_calendar"
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "trial_end_time",
    ignore.case = TRUE
  )

  # LPI min follow-up requires min_followup
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, trt_effect = 0)
      ),
      enrollment_details = list(
        enrollment_distribution = "exponential",
        enrollment_exponential_rate = 1
      ),
      trial_end_details = list(
        type = "last_patient_min_followup"
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "min_followup",
    ignore.case = TRUE
  )

  # event-driven requires a TTE endpoint
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "normal", baseline_mean = 0, sd = 1, trt_effect = 0)
      ),
      trial_end_details = list(
        type = "event_driven",
        event_endpoint = "TTE_1",
        target_events = 5
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "time-to-event|TTE",
    ignore.case = TRUE
  )
})


# --- event-driven trial end -------------------------------------------------

testthat::test_that("event-driven trial end: target events reached and analysis population identified", {
  skip_if_missing_pkgs()

  obj <- makeData(
    correlation_matrix = NULL,
    sample_size_per_group = 500,
    SEED = 7,
    endpoint_details = list(ep_tte_fast),
    enrollment_details = list(
      enrollment_distribution = "exponential",
      enrollment_exponential_rate = 20
    ),
    followup_details = list(
      max_followup = 156
    ),
    trial_end_details = list(
      type = "event_driven",
      event_endpoint = "TTE_1",
      target_events = 150
    ),
    non_fatal_censors_fatal = FALSE,
    target_correlation = FALSE
  )

  d <- obj$data
  tc <- obj$meta$trial_calendar

  testthat::expect_equal(nrow(d), 1000L)

  testthat::expect_true("enrollTime" %in% names(d))
  testthat::expect_true("availableFollowup" %in% names(d))

  testthat::expect_equal(tc$trial_end_reason, "event_target")
  testthat::expect_true(isTRUE(tc$event_target_reached))
  testthat::expect_true(tc$n_events_at_end >= 150)

  expected_fu <- pmin(156, pmax(0, tc$trial_end_time - d$enrollTime))

  testthat::expect_equal(d$availableFollowup, expected_fu, tolerance = 1e-10)
  testthat::expect_true(all(d$TTE_1 <= d$availableFollowup + 1e-12))

  analysis_idx <- d$availableFollowup > 0

  testthat::expect_true(any(analysis_idx))
  testthat::expect_true(sum(analysis_idx) <= nrow(d))

  # Rows after trial end are generated as planned subjects but have no follow-up.
  idx0 <- which(d$availableFollowup == 0)

  if (length(idx0) > 0) {
    testthat::expect_true(all(abs(d$TTE_1[idx0]) < 1e-12))
    testthat::expect_true(all(d$Status_1[idx0] == 0))
    testthat::expect_true(all(d$enrollTime[idx0] >= tc$trial_end_time - 1e-12))
  }

  # Analysis population should be the positive-follow-up population.
  d_analysis <- d[analysis_idx, , drop = FALSE]

  testthat::expect_true(all(d_analysis$availableFollowup > 0))
  testthat::expect_true(sum(d_analysis$Status_1) >= 150)
})


testthat::test_that("event-driven trial end: require_min_followup delays final analysis", {
  skip_if_missing_pkgs()

  min_followup <- 52
  max_followup <- 156

  obj <- makeData(
    correlation_matrix = NULL,
    sample_size_per_group = 500,
    SEED = 8,
    endpoint_details = list(ep_tte_fast),
    enrollment_details = list(
      enrollment_distribution = "exponential",
      enrollment_exponential_rate = 25
    ),
    followup_details = list(
      min_followup = min_followup,
      max_followup = max_followup
    ),
    trial_end_details = list(
      type = "event_driven",
      event_endpoint = 1,
      target_events = 100,
      require_min_followup = TRUE
    ),
    non_fatal_censors_fatal = FALSE,
    target_correlation = FALSE
  )

  d <- obj$data
  tc <- obj$meta$trial_calendar

  lpi_min_time <- max(d$enrollTime) + min_followup

  testthat::expect_true(grepl("plus_min_followup", tc$trial_end_reason))
  testthat::expect_true(tc$trial_end_time >= lpi_min_time - 1e-10)

  expected_fu <- pmin(max_followup, pmax(0, tc$trial_end_time - d$enrollTime))

  testthat::expect_equal(d$availableFollowup, expected_fu, tolerance = 1e-10)
  testthat::expect_true(max(d$availableFollowup) <= max_followup + 1e-12)
})


testthat::test_that("event-driven trial end: target-not-reached fallback works", {
  skip_if_missing_pkgs()

  min_followup <- 52
  max_followup <- 104


  obj <- NULL

  testthat::expect_warning(
    obj <- makeData(
      correlation_matrix = NULL,
      sample_size_per_group = 100,
      SEED = 9,
      endpoint_details = list(ep_tte_calendar),
      enrollment_details = list(
        enrollment_distribution = "exponential",
        enrollment_exponential_rate = 10
      ),
      followup_details = list(
        min_followup = min_followup,
        max_followup = max_followup
      ),
      trial_end_details = list(
        type = "event_driven",
        event_endpoint = "TTE_1",
        target_events = 100000,
        target_not_reached = "last_patient_min_followup"
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "Event target was not reached"
  )



  d <- obj$data
  tc <- obj$meta$trial_calendar

  expected_trial_end <- max(d$enrollTime) + min_followup
  expected_fu <- pmin(max_followup, pmax(0, expected_trial_end - d$enrollTime))

  testthat::expect_false(isTRUE(tc$event_target_reached))
  testthat::expect_equal(
    tc$trial_end_reason,
    "last_patient_min_followup_target_not_reached"
  )
  testthat::expect_equal(tc$trial_end_time, expected_trial_end, tolerance = 1e-10)
  testthat::expect_equal(d$availableFollowup, expected_fu, tolerance = 1e-10)
})


# --- arm_mode = control: no trt column --------------------------------------

testthat::test_that("arm_mode='control': omit trt column and behave like K=1", {
  skip_if_missing_pkgs()

  set.seed(10)

  obj <- makeData(
    correlation_matrix = diag(2),
    sample_size_per_group = 2000,
    SEED = 202,
    arm_mode = "control",
    endpoint_details = list(
      list(endpoint_type = "normal", baseline_mean = 5, sd = 2, trt_effect = NULL),
      list(endpoint_type = "binary", baseline_prob = 0.40, trt_effect = NULL)
    ),
    enrollment_details = list(
      enrollment_distribution = "none"
    ),
    followup_details = list(),
    trial_end_details = list(
      type = "none"
    ),
    non_fatal_censors_fatal = FALSE,
    target_correlation = TRUE,
    calibration_control = cal_ctl_fast
  )

  d <- obj$data

  testthat::expect_false("trt" %in% names(d))
  testthat::expect_true(isTRUE(obj$meta$control_only))
  testthat::expect_equal(obj$meta$n_arms, 1L)

  s <- summary(obj)
  testthat::expect_equal(s$n_arms, 1L)

  testthat::expect_true(abs(mean(d$Cont_1) - 5) < 0.12)
  testthat::expect_true(abs(sd(d$Cont_1) - 2) < 0.12)
  testthat::expect_true(abs(mean(d$Bin_1) - 0.40) < 0.03)
})


# --- single-endpoint mode ---------------------------------------------------

testthat::test_that("single-endpoint mode: correlation_matrix=NULL generates one endpoint", {
  skip_if_missing_pkgs()

  set.seed(11)

  obj <- makeData(
    correlation_matrix = NULL,
    sample_size_per_group = 5000,
    SEED = 303,
    arm_mode = "control",
    endpoint_details = list(
      list(
        endpoint_type = "count",
        baseline_mean = 7,
        size = 4,
        p_zero = 0.1,
        trt_effect = NULL
      )
    ),
    enrollment_details = list(
      enrollment_distribution = "none"
    ),
    followup_details = list(),
    trial_end_details = list(
      type = "none"
    ),
    non_fatal_censors_fatal = FALSE,
    target_correlation = FALSE
  )

  d <- obj$data

  testthat::expect_true(isTRUE(obj$meta$single_endpoint_mode))
  testthat::expect_true("Int_1" %in% names(d))
  testthat::expect_false("trt" %in% names(d))

  mu0 <- 7
  p0  <- 0.1
  mean_target <- (1 - p0) * mu0

  testthat::expect_true(abs(mean(d$Int_1) - mean_target) < 0.35)
  testthat::expect_true(mean(d$Int_1 == 0) > 0.05)
})
