# purpose & overview
#
# 1) "marginals and correlation structure look correct (large n, multi-arm)"
#    Generates a mixed-endpoint, multi-arm trial in a large sample and checks that:
#      - observed (within-arm) Pearson correlation matrices are close to the
#        user-supplied target correlation matrix (off-diagonals, within tolerance),
#      - each endpoint’s marginal distribution matches its inputs (via empirical
#        checks and simple regression-style comparisons).
#
# 2) "semi-competing risks: non-fatal cannot occur after fatal (multi-arm)"
#    Validates the semi-competing risks rule: if a non-fatal endpoint is observed,
#    it must occur on or before the fatal time. Ensures “fatal event censors later endpoints.”
#
# 3) "non_fatal_censors_fatal=TRUE: censoring of non-fatal censors other TTE endpoints (multi-arm)"
#    Validates the optional rule that censoring of a *non-fatal* endpoint can be
#    used to censor other TTE endpoints. We check the time truncation implication
#    but do not require the fatal endpoint to be censored (death-induced censoring
#    of non-fatal should not “back-censor” fatal).
#
# 4) "argument checks: required fields and structural rules enforced"
#    Ensures validation is working: required fields must be present and valid,
#    structural rules (fatal ordering, multi-arm length rules, etc.) are enforced,
#    and mutual exclusivity holds for trt_prob vs trt_effect and trt_count vs trt_effect.
#
# 5) "enrollment: admin censoring uses enrollTime-adjusted follow-up (multi-arm)"
#    Validates staggered entry + administrative censoring: each subject’s max follow-up
#    is (admin_censoring - enrollTime), and TTE times are truncated accordingly.
#
# 6) "enrollment argument dependencies are enforced"
#    Ensures enrollment settings are internally consistent: enrollment != none requires
#    administrative_censoring; exponential requires positive rate; piecewise requires
#    increasing cutpoints, correct number of rates, and positive rates.
#
# 7) "arm_mode='control': omit trt column and behave like K=1"
#    When arm_mode='control', only control group data are generated and the returned
#    data omit the `trt` column. The meta should mark control_only=TRUE.
#
# 8) "single-endpoint mode: correlation_matrix=NULL generates one endpoint"
#    If correlation_matrix=NULL, generator behaves as a single-endpoint sampler.
#    Must error if multiple endpoints are supplied. Meta should mark single_endpoint_mode=TRUE.
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
      c(1,2, 0.20),
      c(1,3, 0.10),
      c(1,4, 0.15),
      c(2,3, 0.25),
      c(2,4, 0.05),
      c(3,4, 0.30)
    )
  )
}

offdiag_max_abs <- function(A, B) {
  idx <- which(row(A) != col(A), arr.ind = TRUE)
  max(abs(A[idx] - B[idx]))
}

cal_ctl_fast <- list(
  n_obs = 5000,
  tol = 0.005,
  maxit = 75,
  rho_max = 0.999,
  ensure_cor_mat = TRUE,
  conv_norm_type = "O"
)

get_cor_arm <- function(s, k) {
  nm <- paste0("arm_", k)
  if (!is.null(s$estimated_correlation_by_arm) && !is.null(s$estimated_correlation_by_arm[[nm]])) {
    return(s$estimated_correlation_by_arm[[nm]])
  }
  NULL
}



# --- 1) large-sample marginals + correlation -------------------------------

testthat::test_that("marginals and correlation structure look correct (large n, multi-arm)", {
  skip_if_missing_pkgs()

  set.seed(1)
  cor_mat <- make_cor4()

  # Force full multi-arm behavior explicitly
  # K=3 arms (control + 2 active)
  n_by_arm <- c(5000, 6500, 5500)

  # heteroskedastic SD by arm for continuous endpoint (length K)
  sd_vec <- c(3.0, 4.0, 2.5)

  obj <- makeData(
    correlation_matrix = cor_mat,
    sample_size_per_group = n_by_arm,
    SEED = 123,
    arm_mode = "full",
    endpoint_details = list(
      list(endpoint_type = "normal", baseline_mean = 10, sd = sd_vec,
           trt_effect = c(-2, -1)),
      list(endpoint_type = "binary", baseline_prob = 0.30,
           trt_prob = c(0.25, 0.20)),
      list(endpoint_type = "count", baseline_mean = 8, size = 5, p_zero = 0,
           trt_count = c(10, 9)),
      list(endpoint_type = "tte", baseline_rate = 1/24, censoring_rate = 0,
           fatal_event = FALSE, trt_effect = c(log(0.80), log(0.90)))
    ),
    enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
    non_fatal_censors_fatal = FALSE,
    target_correlation = TRUE,
    calibration_control = cal_ctl_fast
  )

  s <- summary(obj)
  d <- obj$data

  # --- sample size per arm checks
  testthat::expect_true("trt" %in% names(d))
  tab <- table(d$trt)
  testthat::expect_equal(as.integer(tab["0"]), n_by_arm[1])
  testthat::expect_equal(as.integer(tab["1"]), n_by_arm[2])
  testthat::expect_equal(as.integer(tab["2"]), n_by_arm[3])

  # --- correlations close to target (off-diagonals)
  # allow a bit looser tolerance due to mixed types + heteroskedasticity + calibration stochasticity
  for (k in 0:2) {
    cor_k <- get_cor_arm(s, k)
    testthat::expect_true(!is.null(cor_k))
    testthat::expect_true(offdiag_max_abs(cor_k, cor_mat) < 0.07)
  }

  # --- continuous means
  testthat::expect_true(abs(mean(d$Cont_1[d$trt == 0]) - 10) < 0.25)
  testthat::expect_true(abs(mean(d$Cont_1[d$trt == 1]) -  8) < 0.30)
  testthat::expect_true(abs(mean(d$Cont_1[d$trt == 2]) -  9) < 0.30)

  # new tests to test for exactness (when updating algo cadence)
  testthat::expect_true(s$continuous$est_resid_sd[1] - 3.019113 < 0.0005)
  testthat::expect_true(s$tte$exp_rate[2]- 0.03755223 < 0.0005)


  # --- continuous SDs (heteroskedastic)
  testthat::expect_true(abs(sd(d$Cont_1[d$trt == 0]) - sd_vec[1]) < 0.30)
  testthat::expect_true(abs(sd(d$Cont_1[d$trt == 1]) - sd_vec[2]) < 0.35)
  testthat::expect_true(abs(sd(d$Cont_1[d$trt == 2]) - sd_vec[3]) < 0.30)

  # --- binary probs by arm (trt_prob targets)
  p0 <- mean(d$Bin_1[d$trt == 0] == 1)
  p1 <- mean(d$Bin_1[d$trt == 1] == 1)
  p2 <- mean(d$Bin_1[d$trt == 2] == 1)
  testthat::expect_true(abs(p0 - 0.30) < 0.02)
  testthat::expect_true(abs(p1 - 0.25) < 0.03)
  testthat::expect_true(abs(p2 - 0.20) < 0.03)

  # --- count means by arm (trt_count targets)
  m0 <- mean(d$Int_1[d$trt == 0])
  m1 <- mean(d$Int_1[d$trt == 1])
  m2 <- mean(d$Int_1[d$trt == 2])
  testthat::expect_true(abs(m0 -  8) < 0.40)
  testthat::expect_true(abs(m1 - 10) < 0.55)
  testthat::expect_true(abs(m2 -  9) < 0.55)

  # --- TTE direction check: HR<1 => larger mean time
  t0 <- mean(d$TTE_1[d$trt == 0])
  t1 <- mean(d$TTE_1[d$trt == 1])
  t2 <- mean(d$TTE_1[d$trt == 2])
  testthat::expect_true(t1 > t0)
  testthat::expect_true(t2 > t0)
})



# --- 2) semi-competing risks ------------------------------------------------

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
      list(endpoint_type = "tte", baseline_rate = 1/20, censoring_rate = 1/80,
           fatal_event = TRUE,  trt_effect = c(log(0.9), log(0.95))),
      list(endpoint_type = "tte", baseline_rate = 1/10, censoring_rate = 1/80,
           fatal_event = FALSE, trt_effect = c(0, 0))
    ),
    enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
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
  cor_mat <- corr_make(num_endpoints = 2, values = rbind(c(1,2,0.1)))
  n_by_arm <- c(3000, 3000, 3000)

  # fatal first to satisfy ordering checks
  obj <- makeData(
    correlation_matrix = cor_mat,
    sample_size_per_group = n_by_arm,
    SEED = 999,
    arm_mode = "full",
    endpoint_details = list(
      list(endpoint_type = "tte", baseline_rate = 1/30, censoring_rate = 1/12,
           fatal_event = TRUE,  trt_effect = c(0,0)),
      list(endpoint_type = "tte", baseline_rate = 1/25, censoring_rate = 1/12,
           fatal_event = FALSE, trt_effect = c(0,0))
    ),
    enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
    non_fatal_censors_fatal = TRUE,
    target_correlation = FALSE
  )

  d <- obj$data

  # non-fatal is TTE_2 with indicator Status_2
  idx_cens_nf <- which(d$Status_2 == 0)

  if (length(idx_cens_nf) > 0) {
    # censoring time of non-fatal should truncate the other endpoint time
    testthat::expect_true(all(d$TTE_1[idx_cens_nf] <= d$TTE_2[idx_cens_nf] + 1e-12))
  }
})



# --- 3) argument requirement logic ------------------------------------------

testthat::test_that("argument checks: required fields and structural rules enforced", {

  # single-endpoint mode requires exactly 1 endpoint
  testthat::expect_error(
    makeData(
      correlation_matrix = NULL,
      sample_size_per_group = 10,
      SEED = 1,
      endpoint_details = list(
        list(endpoint_type = "normal", baseline_mean = 0, sd = 1),
        list(endpoint_type = "binary", baseline_prob = 0.3)
      ),
      enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE,
      arm_mode = "control"
    ),
    "exactly one endpoint|single.*endpoint|correlation_matrix",
    ignore.case = TRUE
  )

  cor_mat <- diag(1)

  # continuous missing sd
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = 10,
      SEED = 1,
      endpoint_details = list(list(endpoint_type = "normal", baseline_mean = 0, trt_effect = 0)),
      enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE,
      arm_mode = "full"
    ),
    "requires.*sd|sd",
    ignore.case = TRUE
  )

  # binary baseline prob out of bounds
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = 10,
      SEED = 1,
      endpoint_details = list(list(endpoint_type = "binary", baseline_prob = 1.0, trt_effect = 0)),
      enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE,
      arm_mode = "full"
    ),
    "baseline_prob",
    ignore.case = TRUE
  )

  # count requires size
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = 10,
      SEED = 1,
      endpoint_details = list(list(endpoint_type = "count", baseline_mean = 5, trt_effect = 0)),
      enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE,
      arm_mode = "full"
    ),
    "requires.*size|size",
    ignore.case = TRUE
  )

  # fatal must be first or first+second among TTE endpoints
  cor_mat2 <- diag(3)
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat2,
      sample_size_per_group = c(10,10,10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, fatal_event = FALSE, trt_effect = c(0,0)),
        list(endpoint_type = "tte", baseline_rate = 0.1, fatal_event = FALSE, trt_effect = c(0,0)),
        list(endpoint_type = "tte", baseline_rate = 0.1, fatal_event = TRUE,  trt_effect = c(0,0))
      ),
      enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "Terminal|fatal|first",
    ignore.case = TRUE
  )

  # mutual exclusivity: binary trt_prob and trt_effect cannot both be set
  testthat::expect_error(
    makeData(
      correlation_matrix = diag(1),
      sample_size_per_group = c(10,10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "binary", baseline_prob = 0.3, trt_effect = 0.1, trt_prob = 0.2)
      ),
      enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "trt_prob|trt_effect|only one|exclusive",
    ignore.case = TRUE
  )

  # mutual exclusivity: count trt_count and trt_effect cannot both be set
  testthat::expect_error(
    makeData(
      correlation_matrix = diag(1),
      sample_size_per_group = c(10,10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "count", baseline_mean = 5, size = 2, trt_effect = 0.1, trt_count = 6)
      ),
      enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "trt_count|trt_effect|only one|exclusive",
    ignore.case = TRUE
  )
})


testthat::test_that("argument checks: sd length rules + sample_size length rules enforced", {

  # With arm_mode="full" and trt_effect length 2, K=3
  cor_mat <- diag(1)

  # OK: sd length K (3), sample_size length K (3)
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 12, 11),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "normal", baseline_mean = 0, sd = c(1, 1.2, 0.8),
             trt_effect = c(0.1, 0.2))
      ),
      enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    NA
  )

  # FAIL: sd length not 1 or K
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 12, 11),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "normal", baseline_mean = 0, sd = c(1, 2),
             trt_effect = c(0.1, 0.2))
      ),
      enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "sd.*length 1 or K|sd.*length",
    ignore.case = TRUE
  )

  # FAIL: sample_size_per_group length not 1 or K
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10, 12), # invalid length
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "normal", baseline_mean = 0, sd = 1,
             trt_effect = c(0.1, 0.2))
      ),
      enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "sample_size_per_group|length",
    ignore.case = TRUE
  )
})



# --- 4) enrollment behavior -------------------------------------------------

testthat::test_that("enrollment: admin censoring uses enrollTime-adjusted follow-up (multi-arm)", {
  skip_if_missing_pkgs()

  set.seed(4)
  cor_mat <- diag(1)

  admin <- 12

  obj <- makeData(
    correlation_matrix = cor_mat,
    sample_size_per_group = c(2500, 2500, 2500),
    SEED = 88,
    arm_mode = "full",
    endpoint_details = list(
      list(endpoint_type = "tte", baseline_rate = 1/6, censoring_rate = 0,
           fatal_event = FALSE, trt_effect = c(0,0))
    ),
    enrollment_details = list(
      administrative_censoring = admin,
      enrollment_distribution  = "uniform"
    ),
    non_fatal_censors_fatal = FALSE,
    target_correlation = FALSE
  )

  d <- obj$data

  testthat::expect_true("enrollTime" %in% names(d))
  testthat::expect_true(all(d$enrollTime >= 0))
  testthat::expect_true(all(d$enrollTime <= admin + 1e-12))

  max_fu <- pmax(0, admin - d$enrollTime)
  testthat::expect_true(all(d$TTE_1 <= max_fu + 1e-12))

  idx0 <- which(max_fu == 0)
  if (length(idx0) > 0) {
    testthat::expect_true(all(abs(d$TTE_1[idx0]) < 1e-12))
    testthat::expect_true(all(d$Status_1[idx0] == 0))
  }
})


testthat::test_that("enrollment argument dependencies are enforced", {
  cor_mat <- diag(1)

  # enrollment != none requires administrative_censoring
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10,10,10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, trt_effect = c(0,0))
      ),
      enrollment_details = list(
        administrative_censoring = NULL,
        enrollment_distribution  = "uniform"
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "administrative_censoring",
    ignore.case = TRUE
  )

  # exponential enrollment requires positive rate
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10,10,10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, trt_effect = c(0,0))
      ),
      enrollment_details = list(
        administrative_censoring = 10,
        enrollment_distribution  = "exponential",
        enrollment_exponential_rate = NULL
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "exponential_rate",
    ignore.case = TRUE
  )

  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10,10,10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, trt_effect = c(0,0))
      ),
      enrollment_details = list(
        administrative_censoring = 10,
        enrollment_distribution  = "exponential",
        enrollment_exponential_rate = -1
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "positive",
    ignore.case = TRUE
  )

  # piecewise enrollment checks
  testthat::expect_error(
    makeData(
      correlation_matrix = cor_mat,
      sample_size_per_group = c(10,10,10),
      SEED = 1,
      arm_mode = "full",
      endpoint_details = list(
        list(endpoint_type = "tte", baseline_rate = 0.1, trt_effect = c(0,0))
      ),
      enrollment_details = list(
        administrative_censoring = 10,
        enrollment_distribution  = "piecewise",
        piecewise_enrollment_cutpoints = c(0, 5, 4),
        piecewise_enrollment_rates = c(0.1, 0.2)
      ),
      non_fatal_censors_fatal = FALSE,
      target_correlation = FALSE
    ),
    "increasing",
    ignore.case = TRUE
  )
})



# --- 5) arm_mode = control: no trt column ----------------------------------

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
    enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
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



# --- 6) single-endpoint mode (correlation_matrix=NULL) ----------------------

testthat::test_that("single-endpoint mode: correlation_matrix=NULL generates one endpoint", {
  skip_if_missing_pkgs()

  set.seed(11)

  obj <- makeData(
    correlation_matrix = NULL,
    sample_size_per_group = 5000,
    SEED = 303,
    arm_mode = "control",
    endpoint_details = list(
      list(endpoint_type = "count", baseline_mean = 7, size = 4, p_zero = 0.1, trt_effect = NULL)
    ),
    enrollment_details = list(administrative_censoring = NULL, enrollment_distribution = "none"),
    non_fatal_censors_fatal = FALSE,
    target_correlation = FALSE
  )

  d <- obj$data

  testthat::expect_true(isTRUE(obj$meta$single_endpoint_mode))
  testthat::expect_true("Int_1" %in% names(d))
  testthat::expect_false("trt" %in% names(d))

  # Under the current qzinb_mixture() implementation:
  # E[Y] = (1 - p0) * mu
  mu0 <- 7
  p0  <- 0.1
  mean_target <- (1 - p0) * mu0  # 6.3

  testthat::expect_true(abs(mean(d$Int_1) - mean_target) < 0.35)

  # sanity check: some zeros should exist
  testthat::expect_true(mean(d$Int_1 == 0) > 0.05)
})

