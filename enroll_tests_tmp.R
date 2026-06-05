# ==============================================================================
# Manual tests for updated makeData() calendar features
# ==============================================================================

library(endpoints)

# ------------------------------------------------------------------------------
# Small assertion helpers
# ------------------------------------------------------------------------------

check_true <- function(x, msg) {
  if (!isTRUE(x)) stop("FAILED: ", msg, call. = FALSE)
  message("PASSED: ", msg)
}

check_equal <- function(x, y, msg, tol = 1e-8) {
  ok <- isTRUE(all.equal(x, y, tolerance = tol))
  if (!ok) stop("FAILED: ", msg, "\n  Got: ", x, "\n  Expected: ", y, call. = FALSE)
  message("PASSED: ", msg)
}

# ------------------------------------------------------------------------------
# Endpoint definitions
# ------------------------------------------------------------------------------

ep_cont <- list(
  endpoint_type = "continuous",
  baseline_mean = 10,
  sd            = 2,
  trt_effect    = -1
)

ep_tte <- list(
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

# ==============================================================================
# 1. Baseline: no enrollment, no follow-up details, no trial end
# ==============================================================================

message("\n--- Test 1: baseline makeData() behavior ---")

sim1 <- makeData(
  correlation_matrix    = NULL,
  sample_size_per_group = 100,
  SEED                  = 1,
  endpoint_details      = list(ep_tte)
)

dat1 <- as.data.frame(sim1)

check_true("TTE_1" %in% names(dat1), "TTE_1 exists")
check_true("Status_1" %in% names(dat1), "Status_1 exists")
check_true(!"enrollTime" %in% names(dat1), "enrollTime absent by default")
check_true(!"availableFollowup" %in% names(dat1), "availableFollowup absent by default")

print(sim1)
summary(sim1)

# ==============================================================================
# 2. Exponential enrollment = homogeneous Poisson-process accrual
# ==============================================================================

message("\n--- Test 2: exponential / Poisson enrollment ---")

sim2 <- makeData(
  correlation_matrix    = NULL,
  sample_size_per_group = 100,
  SEED                  = 2,
  endpoint_details      = list(ep_cont),
  enrollment_details    = list(
    enrollment_distribution     = "exponential",
    enrollment_exponential_rate = 10
  )
)

dat2 <- as.data.frame(sim2)

check_true("enrollTime" %in% names(dat2), "enrollTime exists with stochastic enrollment")
check_true(min(dat2$enrollTime) >= 0, "enrollTime non-negative")
check_true(max(dat2$enrollTime) > 0, "enrollTime has positive range")

summary(dat2$enrollTime)
hist(dat2$enrollTime, main = "Poisson-process enrollment times", xlab = "Enrollment time")

# ==============================================================================
# 3. max_followup only
# ==============================================================================

message("\n--- Test 3: max_followup only ---")

sim3 <- makeData(
  correlation_matrix    = NULL,
  sample_size_per_group = 200,
  SEED                  = 3,
  endpoint_details      = list(ep_tte),
  followup_details      = list(
    max_followup = 52
  )
)

dat3 <- as.data.frame(sim3)

check_true("availableFollowup" %in% names(dat3), "availableFollowup exists with max_followup")
check_true(all(dat3$availableFollowup == 52), "availableFollowup equals max_followup")
check_true(max(dat3$TTE_1) <= 52 + 1e-8, "TTE capped by max_followup")

table(dat3$Status_1)
summary(dat3)

# ==============================================================================
# 4. last_patient_min_followup
# ==============================================================================

message("\n--- Test 4: last_patient_min_followup ---")

sim4 <- makeData(
  correlation_matrix    = NULL,
  sample_size_per_group = 150,
  SEED                  = 4,
  endpoint_details      = list(ep_tte),
  enrollment_details    = list(
    enrollment_distribution     = "exponential",
    enrollment_exponential_rate = 8
  ),
  followup_details = list(
    min_followup = 52,
    max_followup = 156
  ),
  trial_end_details = list(
    type = "last_patient_min_followup"
  )
)

dat4 <- as.data.frame(sim4)
tc4  <- sim4$meta$trial_calendar

check_true("enrollTime" %in% names(dat4), "enrollTime exists")
check_true("availableFollowup" %in% names(dat4), "availableFollowup exists")
check_equal(tc4$trial_end_reason, "last_patient_min_followup", "trial end reason is LPI min follow-up")

expected_end4 <- max(dat4$enrollTime) + 52

check_equal(tc4$trial_end_time, expected_end4, "trial_end_time = max(enrollTime) + min_followup")
check_true(max(dat4$availableFollowup) <= 156 + 1e-8, "availableFollowup capped by max_followup")
check_true(min(dat4$availableFollowup) <= 52 + 1e-8, "latest enroller has approximately min follow-up")

summary(dat4$enrollTime)
summary(dat4$availableFollowup)
print(sim4)

# ==============================================================================
# 5. Fixed calendar trial end
# ==============================================================================

message("\n--- Test 5: fixed_calendar trial end ---")

sim5 <- makeData(
  correlation_matrix    = NULL,
  sample_size_per_group = 150,
  SEED                  = 5,
  endpoint_details      = list(ep_tte),
  enrollment_details    = list(
    enrollment_distribution     = "exponential",
    enrollment_exponential_rate = 5
  ),
  #followup_details = list(
  #  max_followup = 104
  #),
  trial_end_details = list(
    type = "fixed_calendar",
    trial_end_time = 80
  )
)

dat5 <- as.data.frame(sim5)
tc5  <- sim5$meta$trial_calendar

check_equal(tc5$trial_end_reason, "fixed_calendar", "trial end reason fixed_calendar")
check_equal(tc5$trial_end_time, 80, "trial_end_time equals fixed calendar time")

expected_fu5 <- pmin(104, pmax(0, 80 - dat5$enrollTime))

check_true(all(abs(dat5$availableFollowup - expected_fu5) < 1e-8),
           "availableFollowup equals min(max_followup, trial_end_time - enrollTime)")

check_true(max(dat5$TTE_1) <= max(dat5$availableFollowup) + 1e-8,
           "TTE times respect available follow-up")

summary(dat5$availableFollowup)
message("Subjects with zero available follow-up: ", sum(dat5$availableFollowup == 0))

# ==============================================================================
# 6. Piecewise Poisson enrollment
# ==============================================================================

message("\n--- Test 6: piecewise Poisson enrollment ---")

sim6 <- makeData(
  correlation_matrix    = NULL,
  sample_size_per_group = 250,
  SEED                  = 6,
  endpoint_details      = list(ep_tte),
  enrollment_details    = list(
    enrollment_distribution        = "piecewise",
    piecewise_enrollment_cutpoints = c(0, 12, 24, 52),
    piecewise_enrollment_rates     = c(2, 8, 12)
  ),
  followup_details = list(
    min_followup = 52,
    max_followup = 156
  ),
  trial_end_details = list(
    type = "last_patient_min_followup"
  )
)

dat6 <- as.data.frame(sim6)

check_true("enrollTime" %in% names(dat6), "piecewise enrollTime exists")
check_true(max(dat6$enrollTime) > 0, "piecewise enrollment has positive accrual duration")
check_true(max(dat6$availableFollowup) <= 156 + 1e-8, "max follow-up cap respected")

summary(dat6$enrollTime)
summary(dat6$availableFollowup)

table(cut(
  dat6$enrollTime,
  breaks = c(0, 12, 24, 52, Inf),
  include.lowest = TRUE
))

hist(dat6$enrollTime, main = "Piecewise Poisson enrollment times", xlab = "Enrollment time")

# ==============================================================================
# 7. Event-driven trial end
# ==============================================================================

message("\n--- Test 7: event-driven trial end ---")

sim7 <- makeData(
  correlation_matrix    = NULL,
  sample_size_per_group = 500,
  SEED                  = 7,
  endpoint_details      = list(ep_tte_fast),
  enrollment_details    = list(
    enrollment_distribution     = "exponential",
    enrollment_exponential_rate = 20
  ),
  followup_details = list(
    max_followup = 156
  ),
  trial_end_details = list(
    type = "event_driven",
    event_endpoint = "TTE_1",
    target_events = 150
  )
)

dat7 <- as.data.frame(sim7)
tc7  <- sim7$meta$trial_calendar

check_equal(tc7$trial_end_reason, "event_target", "event-driven trial ended by event target")
check_true(isTRUE(tc7$event_target_reached), "event target reached")
check_true(tc7$n_events_at_end >= 150, "number of events at end is at least target")
check_true("availableFollowup" %in% names(dat7), "availableFollowup exists")

print(tc7)
table(dat7$Status_1)
summary(dat7$availableFollowup)

# ==============================================================================
# 8. Event-driven plus minimum follow-up
# ==============================================================================

message("\n--- Test 8: event-driven plus minimum follow-up ---")

sim8 <- makeData(
  correlation_matrix    = NULL,
  sample_size_per_group = 500,
  SEED                  = 8,
  endpoint_details      = list(ep_tte_fast),
  enrollment_details    = list(
    enrollment_distribution     = "exponential",
    enrollment_exponential_rate = 25
  ),
  followup_details = list(
    min_followup = 52,
    max_followup = 156
  ),
  trial_end_details = list(
    type = "event_driven",
    event_endpoint = 1,
    target_events = 100,
    require_min_followup = TRUE
  )
)

dat8 <- as.data.frame(sim8)
tc8  <- sim8$meta$trial_calendar

lpi_min_time8 <- max(dat8$enrollTime) + 52

check_true(grepl("plus_min_followup", tc8$trial_end_reason),
           "trial end reason includes plus_min_followup")

check_true(tc8$trial_end_time >= lpi_min_time8 - 1e-8,
           "trial end occurs no earlier than LPI + min_followup")

print(tc8)
summary(dat8$availableFollowup)

# ==============================================================================
# 9. Event target not reached fallback
# ==============================================================================

message("\n--- Test 9: event target not reached fallback ---")

sim9 <- makeData(
  correlation_matrix    = NULL,
  sample_size_per_group = 100,
  SEED                  = 9,
  endpoint_details      = list(ep_tte),
  enrollment_details    = list(
    enrollment_distribution     = "exponential",
    enrollment_exponential_rate = 10
  ),
  followup_details = list(
    min_followup = 52,
    max_followup = 104
  ),
  trial_end_details = list(
    type = "event_driven",
    event_endpoint = "TTE_1",
    target_events = 100000,
    target_not_reached = "last_patient_min_followup"
  )
)

dat9 <- as.data.frame(sim9)
tc9  <- sim9$meta$trial_calendar

check_true(!isTRUE(tc9$event_target_reached), "event target not reached")
check_equal(tc9$trial_end_reason,
            "last_patient_min_followup_target_not_reached",
            "fallback reason is LPI min-follow-up")

expected_end9 <- max(dat9$enrollTime) + 52

check_equal(tc9$trial_end_time, expected_end9, "fallback trial end time correct")

print(tc9)

# ==============================================================================
# 10. Correlated continuous + TTE endpoint with calendar features
# ==============================================================================

message("\n--- Test 10: correlated continuous + TTE with calendar features ---")

R2 <- corr_make(
  num_endpoints = 2,
  values = rbind(c(1, 2, 0.25))
)

sim10 <- makeData(
  correlation_matrix    = R2,
  sample_size_per_group = 300,
  SEED                  = 10,
  endpoint_details      = list(ep_cont, ep_tte),
  enrollment_details    = list(
    enrollment_distribution     = "exponential",
    enrollment_exponential_rate = 12
  ),
  followup_details = list(
    min_followup = 52,
    max_followup = 156
  ),
  trial_end_details = list(
    type = "last_patient_min_followup"
  ),
  target_correlation = TRUE,
  calibration_control = list(
    n_mc = 5000,
    tol = 0.005,
    maxit = 50,
    rho_cap = 0.999,
    ensure_pd = TRUE,
    conv_norm_type = "F"
  )
)

dat10 <- as.data.frame(sim10)

check_true(all(c("Cont_1", "TTE_1", "Status_1", "enrollTime", "availableFollowup") %in% names(dat10)),
           "expected columns exist in correlated calendar simulation")

summary(sim10)

# ==============================================================================
# 11. Expected errors
# ==============================================================================

message("\n--- Test 11A: expected error, fixed_calendar without trial_end_time ---")

err1 <- try(
  makeData(
    correlation_matrix    = NULL,
    sample_size_per_group = 50,
    SEED                  = 11,
    endpoint_details      = list(ep_tte),
    trial_end_details = list(
      type = "fixed_calendar"
    )
  ),
  silent = TRUE
)

check_true(inherits(err1, "try-error"), "fixed_calendar without trial_end_time errors")


message("\n--- Test 11B: expected error, LPI rule without min_followup ---")

err2 <- try(
  makeData(
    correlation_matrix    = NULL,
    sample_size_per_group = 50,
    SEED                  = 12,
    endpoint_details      = list(ep_tte),
    enrollment_details = list(
      enrollment_distribution     = "exponential",
      enrollment_exponential_rate = 10
    ),
    trial_end_details = list(
      type = "last_patient_min_followup"
    )
  ),
  silent = TRUE
)

check_true(inherits(err2, "try-error"), "last_patient_min_followup without min_followup errors")


message("\n--- Test 11C: expected error, event-driven without TTE endpoint ---")

err3 <- try(
  makeData(
    correlation_matrix    = NULL,
    sample_size_per_group = 50,
    SEED                  = 13,
    endpoint_details      = list(ep_cont),
    trial_end_details = list(
      type = "event_driven",
      event_endpoint = "TTE_1",
      target_events = 10
    )
  ),
  silent = TRUE
)

check_true(inherits(err3, "try-error"), "event_driven without TTE endpoint errors")

# ==============================================================================
# Done
# ==============================================================================

message("\nAll new calendar/enrollment feature tests completed.")
