# Summarize simulated data from a `makeDataSim` object

`summary.makeDataSim()` computes arm-specific diagnostic summaries for a
`"makeDataSim"` object returned by
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md).

## Usage

``` r
# S3 method for class 'makeDataSim'
summary(object, ...)
```

## Arguments

- object:

  An object of class `"makeDataSim"`, typically created by
  [`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md).

- ...:

  Currently unused. Included for S3 method compatibility.

## Value

A list of class `"summary.makeDataSim"` with components:

- target_correlation:

  The target correlation matrix supplied to
  [`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md),
  or `NULL` in single-endpoint mode.

- estimated_correlation_by_arm:

  A named list of empirical arm-specific correlation matrices.

- continuous:

  A data frame of continuous-endpoint summaries, or `NULL` if no
  continuous endpoints were simulated.

- binary:

  A data frame of binary-endpoint summaries, or `NULL` if no binary
  endpoints were simulated.

- count:

  A data frame of count-endpoint summaries, or `NULL` if no count
  endpoints were simulated.

- tte:

  A data frame of time-to-event summaries, or `NULL` if no TTE endpoints
  were simulated.

- n_arms:

  The total number of treatment arms represented in the simulated
  dataset.

## Details

The summary includes:

- the target correlation matrix supplied to
  [`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md),

- empirical endpoint correlations within each treatment arm, and

- endpoint-specific marginal summaries for continuous, binary, count,
  and time-to-event outcomes.

These summaries are intended as a quick validation tool for checking
that the simulated dataset is broadly consistent with the requested
data-generating parameters.

## General behavior

The summary method extracts the simulated dataset and metadata stored in
the `"makeDataSim"` object, then computes:

1.  the requested target correlation matrix, as stored in
    `object$meta$correlation_matrix`,

2.  the observed Pearson correlation matrix among endpoint columns
    within each treatment arm, and

3.  endpoint-type-specific summaries based on simple fitted models or
    direct descriptive estimators.

## Arm-specific empirical correlation

For each study arm, the method computes the Pearson correlation matrix
of the simulated endpoint columns listed in
`object$meta$endpoint_names`. These are returned as a named list in
`$estimated_correlation_by_arm`, with elements named `"arm_0"`,
`"arm_1"`, and so on.

Correlations are computed using
[`cor()`](https://rdrr.io/r/stats/cor.html) and rounded to 3 decimal
places.

## Continuous endpoints

For each continuous endpoint:

- a linear model is fit, `lm(y ~ trt)` in multi-arm settings and
  `lm(y ~ 1)` otherwise,

- the control-group intercept is used as the estimated baseline mean,

- treatment coefficients are reported as estimated mean shifts, and

- arm-specific residual SDs are computed from the residuals of the
  fitted model.

The returned table includes:

- endpoint:

  Endpoint name, for example `Cont_1`.

- arm:

  Arm index.

- input_baseline_mean:

  Encoded control-group mean.

- input_sd:

  Encoded SD for that arm.

- input_trt_effect:

  Encoded treatment effect for that arm.

- est_baseline_mean:

  Estimated control-group mean from the fitted model.

- est_trt_effect:

  Estimated mean shift for that arm versus control.

- est_resid_sd:

  Residual SD within that arm.

## Binary endpoints

For each binary endpoint:

- a logistic regression model is fit using
  `glm(..., family = binomial())`,

- the control-group intercept is interpreted on the logit scale,

- treatment coefficients are reported as estimated log-odds ratios, and

- observed arm-specific event probabilities are computed directly as
  sample means.

The returned table includes:

- endpoint:

  Endpoint name, for example `Bin_1`.

- arm:

  Arm index.

- input_baseline_prob:

  Encoded control-group probability.

- input_trt_logOR:

  Encoded treatment effect on the log-odds scale.

- input_trt_prob:

  Encoded arm-specific probability.

- est_baseline_prob:

  Estimated control-group probability from the fitted model.

- est_trt_logOR:

  Estimated treatment log-odds ratio for that arm versus control.

- est_prob:

  Observed arm-specific event proportion.

## Count endpoints

For each count endpoint:

- a negative binomial model is fit using
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html),

- the control-group intercept is exponentiated to obtain the estimated
  baseline mean,

- treatment coefficients are reported as estimated log rate-ratios,

- the fitted dispersion parameter is reported, and

- observed arm-specific means and structural zero proportions are
  computed directly.

The returned table includes:

- endpoint:

  Endpoint name, for example `Int_1`.

- arm:

  Arm index.

- input_baseline_mean:

  Specified control-group mean.

- input_trt_logRR:

  Specified treatment effect on the log rate-ratio scale.

- input_trt_mean:

  Specified arm-specific mean count.

- input_size:

  Specified negative-binomial size parameter.

- input_p_zero:

  Specified structural zero probability.

- est_baseline_mean:

  Estimated control-group mean from the fitted model.

- est_trt_logRR:

  Estimated treatment log rate-ratio for that arm versus control.

- est_size:

  Estimated negative-binomial size parameter.

- obs_mean:

  Observed arm-specific mean count.

- obs_p0:

  Observed proportion of zeros in that arm.

## Time-to-event endpoints

For each time-to-event endpoint:

- a Cox proportional hazards model is fit using
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html),

- treatment coefficients are reported as estimated log hazard-ratios,

- observed arm-specific event proportions are reported, and

- an arm-specific exponential maximum likelihood estimate of the event
  rate is computed as \\\hat\lambda = \sum_i \delta_i / \sum_i t_i\\,
  where \\t_i\\ is observed follow-up time and \\\delta_i\\ is the event
  indicator.

The returned table includes:

- endpoint:

  Endpoint name, for example `TTE_1`.

- arm:

  Arm index.

- censor_col:

  Corresponding event-indicator column name, for example `Status_1`.

- input_baseline_rate:

  Requested control-group exponential event rate.

- input_trt_logHR:

  Requested treatment effect on the log hazard-ratio scale.

- est_trt_logHR:

  Estimated log hazard-ratio from the Cox model.

- input_trt_HR:

  Requested treatment effect on the hazard-ratio scale.

- est_trt_HR:

  Estimated hazard-ratio from the Cox model.

- obs_event_rate:

  Observed event proportion in that arm.

- exp_rate:

  Arm-specific exponential MLE event-rate estimate.

## Interpretation

Because simulations are finite-sample and may involve nonlinear marginal
transformations, censoring, administrative censoring, and
fatal/non-fatal event logic, the estimated summaries will generally not
match the requested inputs exactly. The summary output is therefore best
viewed as a diagnostic check rather than an exact recovery target.

## See also

[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
for the main simulation routine.

`print.summary.makeDataSim` for the print method.

[`plot.makeDataSim`](https://boehringer-ingelheim.github.io/endpoints/reference/plot.makeDataSim.md)
for graphical diagnostics.

[`corr_make`](https://boehringer-ingelheim.github.io/endpoints/reference/corr_make.md)
for building correlation matrices.

## Examples

``` r
## Continuous + binary example
ep1 <- list(
  endpoint_type = "continuous",
  baseline_mean = 10,
  sd            = 2,
  trt_effect    = -1
)

ep2 <- list(
  endpoint_type = "binary",
  baseline_prob = 0.30,
  trt_prob      = 0.45
)

R2 <- corr_make(
  num_endpoints = 2,
  values = rbind(c(1, 2, 0.2))
)

sim_obj <- makeData(
  correlation_matrix    = R2,
  sample_size_per_group = 500,
  SEED                  = 1,
  endpoint_details      = list(ep1, ep2)
)

ss <- summary(sim_obj)

## Display top-level structure
ss$n_arms
#> [1] 2
ss$target_correlation
#>      [,1] [,2]
#> [1,]  1.0  0.2
#> [2,]  0.2  1.0
ss$estimated_correlation_by_arm
#> $arm_0
#>        Cont_1 Bin_1
#> Cont_1  1.000 0.149
#> Bin_1   0.149 1.000
#> 
#> $arm_1
#>        Cont_1 Bin_1
#> Cont_1  1.000 0.143
#> Bin_1   0.143 1.000
#> 

## Endpoint-specific summaries
ss$continuous
#>   endpoint arm input_baseline_mean input_sd input_trt_effect est_baseline_mean
#> 1   Cont_1   0                  10        2                0           10.2236
#> 2   Cont_1   1                  10        2               -1           10.2236
#>   est_trt_effect est_resid_sd
#> 1       0.000000     1.945490
#> 2      -1.341728     1.975195
ss$binary
#>   endpoint arm input_baseline_prob input_trt_logOR input_trt_prob
#> 1    Bin_1   0                 0.3       0.0000000           0.30
#> 2    Bin_1   1                 0.3       0.6466272           0.45
#>   est_baseline_prob est_trt_logOR est_prob
#> 1             0.332     0.0000000    0.332
#> 2             0.332     0.4823075    0.446

## Time-to-event example
ep_tte <- list(
  endpoint_type  = "tte",
  baseline_rate  = 1 / 24,
  trt_effect     = log(0.8),
  censoring_rate = 1 / 216,
  fatal_event    = TRUE
)

sim_tte <- makeData(
  correlation_matrix    = NULL,
  sample_size_per_group = 1000,
  SEED                  = 2,
  endpoint_details      = list(ep_tte)
)

summary(sim_tte)$tte
#>   endpoint arm censor_col input_baseline_rate input_trt_logHR input_trt_HR
#> 1    TTE_1   0   Status_1          0.04166667       0.0000000          1.0
#> 2    TTE_1   1   Status_1          0.04166667      -0.2231436          0.8
#>   est_trt_logHR est_trt_HR obs_event_rate   exp_rate
#> 1     0.0000000  1.0000000          0.884 0.04030938
#> 2    -0.1968895  0.8212814          0.871 0.03303419
```
