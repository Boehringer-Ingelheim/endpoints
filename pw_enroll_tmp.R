#library(endpoints)
#library(dplyr)
#library(ggplot2)

# Continuous endpoint just so makeData() has something simple to simulate.
# Because trt_effect is supplied, this will generate a two-arm trial.
ep_cont <- list(
  endpoint_type = "continuous",
  baseline_mean = 10,
  sd = 2,
  trt_effect = -1
)

# Piecewise enrollment example:
#   0-12:  1 subject per time unit
#   12-24: 5 subjects per time unit
#   24-52: 10 subjects per time unit
#
# With sample_size_per_group = 300 and two arms, total N = 600.
sim_piece <- makeData(
  correlation_matrix = NULL,
  sample_size_per_group = 300,
  SEED = 123,
  endpoint_details = list(ep_cont),
  enrollment_details = list(
    enrollment_distribution        = "piecewise",
    piecewise_enrollment_cutpoints = c(0, 12, 24, 52),
    piecewise_enrollment_rates     = c(1, 5, 10)
  )
)

dat_piece <- as.data.frame(sim_piece)

# Basic checks
nrow(dat_piece)
summary(dat_piece$enrollTime)

dat_piece %>%
  summarise(
    n = n(),
    min_enroll = min(enrollTime),
    q1_enroll = quantile(enrollTime, 0.25),
    median_enroll = median(enrollTime),
    mean_enroll = mean(enrollTime),
    q3_enroll = quantile(enrollTime, 0.75),
    max_enroll = max(enrollTime)
  )


ggplot(dat_piece, aes(x = enrollTime)) +
  geom_histogram(binwidth = 2, boundary = 0, color = "white") +
  geom_vline(xintercept = c(12, 24, 52), linetype = "dashed") +
  labs(
    title = "Piecewise enrollment times",
    x = "Enrollment time",
    y = "Number enrolled"
  )



enroll_curve <- dat_piece %>%
  arrange(enrollTime) %>%
  mutate(
    cumulative_n = row_number()
  )

ggplot(enroll_curve, aes(x = enrollTime, y = cumulative_n)) +
  geom_step() +
  geom_vline(xintercept = c(12, 24, 52), linetype = "dashed") +
  labs(
    title = "Cumulative enrollment curve",
    x = "Enrollment time")



cutpoints <- c(0, 12, 24, 52)
rates <- c(1, 5, 10)

interval_counts <- dat_piece %>%
  mutate(
    enroll_interval = cut(
      enrollTime,
      breaks = c(cutpoints, Inf),
      right = FALSE,
      include.lowest = TRUE,
      labels = c("[0,12)", "[12,24)", "[24,52)", "[52,Inf)")
    )
  ) %>%
  count(enroll_interval, name = "n_enrolled")

interval_counts



interval_summary <- interval_counts %>%
  mutate(
    interval_length = c(diff(cutpoints), NA_real_),
    target_rate = c(rates, tail(rates, 1)),
    expected_n = interval_length * target_rate,
    observed_rate = n_enrolled / interval_length
  )

interval_summary


