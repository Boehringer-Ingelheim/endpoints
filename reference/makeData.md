# Simulate trial data with mixed endpoint types, optional Gaussian-copula dependence, and enrollment/censoring features

`makeData()` simulates subject-level trial data with one or more
endpoints. Supported endpoint types include continuous, binary, count,
and time-to-event. When multiple endpoints are supplied, dependence may
be induced through a Gaussian copula (NORTA-style construction), with
optional numerical calibration to approximately match a target Pearson
correlation matrix on the observed endpoint scale.

## Usage

``` r
makeData(
  correlation_matrix,
  SEED = NULL,
  sample_size_per_group,
  endpoint_details,
  enrollment_details = list(),
  non_fatal_censors_fatal = FALSE,
  target_correlation = TRUE,
  arm_mode = c("auto", "full", "control"),
  calibration_control = list(n_mc = 10000, tol = 0.001, maxit = 100, rho_cap = 0.999,
    ensure_pd = TRUE, conv_norm_type = "F")
)
```

## Arguments

- correlation_matrix:

  A numeric correlation matrix specifying the dependence structure
  across endpoints (Pearson correlation on the *latent Gaussian* scale).
  If `NULL`, `makeData()` enters *single-endpoint mode*, and exactly one
  endpoint must be supplied in `endpoint_details`. Typically created
  with
  [`corr_make`](https://boehringer-ingelheim.github.io/endpoints/reference/corr_make.md).

- SEED:

  Optional numeric scalar seed used to initialize the random-number
  generator via [`set.seed()`](https://rdrr.io/r/base/Random.html). If
  `NULL`, the current RNG state is used.

- sample_size_per_group:

  Integer scalar or integer vector giving the sample size per arm. If a
  scalar, the same sample size is used for all arms. If a vector, it
  must be of length \\G\\, where \\K\\ is the total number of
  arms(including control).

- endpoint_details:

  A non-empty list of endpoint specification lists, one per endpoint.

  Each sub-list defines the marginal distribution and arm-specific
  parameters for one endpoint. See
  [`endpoint_details`](https://boehringer-ingelheim.github.io/endpoints/reference/endpoint_details.md)
  for the complete schema, supported endpoint types, and detailed
  examples.

- enrollment_details:

  A named list controlling enrollment and administrative follow-up.
  Fields include:

  list("administrative_censoring")

  :   Numeric scalar. If non-`NULL`, follow-up is administratively
      censored at this time.

  list("enrollment_distribution")

  :   Character. One of `"none"`, `"uniform"`, `"exponential"`,
      `"piecewise"`.

  list("enrollment_exponential_rate")

  :   Numeric scalar rate for exponential enrollment if
      `enrollment_distribution="exponential"`.

  list("piecewise_enrollment_cutpoints")

  :   Numeric vector of cutpoints defining intervals for piecewise
      enrollment.

  list("piecewise_enrollment_rates")

  :   Numeric vector of rates (one per interval) for piecewise
      exponential enrollment.

  See
  [`enrollment_details`](https://boehringer-ingelheim.github.io/endpoints/reference/enrollment_details.md)
  for full details and recommended parameterization.

- non_fatal_censors_fatal:

  Logical. Controls semi-competing risks behavior when multiple TTE
  endpoints are present. If `TRUE`, censoring of a non-fatal endpoint
  can censor subsequent TTE endpoints according to the package’s
  semi-competing-risk rules. If `FALSE`, non-fatal censoring does not
  censor other TTE endpoints.

- target_correlation:

  Logical scalar.

  If `TRUE` and `correlation_matrix` is not `NULL`, a numerical
  calibration step is used to search for a latent Gaussian correlation
  matrix whose transformed endpoints approximately match the requested
  observed Pearson correlation matrix.

  If `FALSE`, the supplied `correlation_matrix` is used directly as the
  latent Gaussian correlation matrix.

  This should value should generally be set to `TRUE`, unless there is a
  compelling reason otherwise.

- arm_mode:

  Character string controlling how the number of treatment arms is
  determined. Must be one of:

  list("\\auto\\")

  :   Infer the number of arms from the endpoint specifications. If no
      treatment-specific quantities are supplied, control-only mode is
      used.

  list("\\full\\")

  :   Force a full multi-arm interpretation, even if some endpoint
      specifications omit treatment-specific effects.

  list("\\control\\")

  :   Force control-only mode. In this case, treatment-specific
      arguments such as `trt_effect`, `trt_prob`, and `trt_count` are
      not allowed.

  This should value should generally be set to `auto`.

- calibration_control:

  Named list of controls for the correlation calibration routine used
  when `target_correlation=TRUE`. Common controls include Monte Carlo
  size, tolerance, iteration caps, and constraints to ensure a valid
  correlation matrix. See
  [`calibration_control`](https://boehringer-ingelheim.github.io/endpoints/reference/calibration_control.md)
  for more detail.

## Value

An object of class `"makeDataSim"`.

Internally, this is a list with at least:

- list("data"):

  A `data.frame` containing the simulated dataset.

- list("meta"):

  A metadata list storing endpoint types, endpoint names, arm counts,
  and input settings.

The `data` component includes some subset of:

- list("trt"):

  Treatment arm indicator (omitted in control-only mode).

- list("Cont_1"):

  Continuous endpoints.

- , :

  Continuous endpoints.

- list("Cont_2"):

  Continuous endpoints.

- , :

  Continuous endpoints.

- list():

  Continuous endpoints.

- list("Bin_1"):

  Binary endpoints.

- , :

  Binary endpoints.

- list("Bin_2"):

  Binary endpoints.

- , :

  Binary endpoints.

- list():

  Binary endpoints.

- list("Int_1"):

  Count endpoints.

- , :

  Count endpoints.

- list("Int_2"):

  Count endpoints.

- , :

  Count endpoints.

- list():

  Count endpoints.

- list("TTE_1"):

  Observed time-to-event variables.

- , :

  Observed time-to-event variables.

- list("TTE_2"):

  Observed time-to-event variables.

- , :

  Observed time-to-event variables.

- list():

  Observed time-to-event variables.

- list("Status_1"):

  TTE event indicators (`1 = event`, `0 = censored`).

- , :

  TTE event indicators (`1 = event`, `0 = censored`).

- list("Status_2"):

  TTE event indicators (`1 = event`, `0 = censored`).

- , :

  TTE event indicators (`1 = event`, `0 = censored`).

- list():

  TTE event indicators (`1 = event`, `0 = censored`).

- list("enrollTime"):

  Enrollment times, when stochastic enrollment is used.

## Details

The function also supports:

- multiple treatment arms,

- independent censoring for time-to-event outcomes,

- fatal and non-fatal time-to-event logic (including semi-competing
  risks),

- administrative censoring,

- stochastic enrollment (uniform, exponential, or piecewise
  exponential),

- and can be used to generate longitudinal data.

## Overview

When multiple endpoints are simulated, `makeData()` uses a Gaussian
copula construction. If \\\mathbf{Z} \sim N(\mathbf{0}, \Sigma_Z)\\ is a
latent multivariate normal vector, then each endpoint is generated as
\$\$ X_j = F_j^{-1}\\\Phi(Z_j)\\, \$\$ where \\\Phi\\ is the standard
normal CDF and \\F_j^{-1}\\ is the endpoint-specific quantile function.

For continuous endpoints this yields Gaussian margins; for binary,
count, and time-to-event endpoints, the corresponding marginal quantile
functions are applied to the copula uniforms.

## Treatment arms

The total number of study arms is generally determined from the lengths
of treatment-specific inputs in `endpoint_details` (e.g. the length of
`trt_effect`). %

- `trt_effect` for continuous, binary, count, and TTE endpoints, %

- `trt_prob` for binary endpoints, %

- `trt_count` for count endpoints. %

When treatment arms are present, the output includes a `trt` column
coded as `0, 1, 2, ...`. This can also be controlled via `arm_mode`.

## Administrative censoring and enrollment

If `administrative_censoring` is supplied, all TTE outcomes are
truncated at the maximum available follow-up. If stochastic enrollment
is enabled, each subject receives an `enrollTime`, and maximum
observable follow-up is reduced to \\\mathcal{A} - T_E\\, where
\\\mathcal{A}\\ is the administrative censoring time and \\T_E\\ is the
enrollment time.

## Output object

The returned object has class `"makeDataSim"` and is intended to be used
with:

- [`print()`](https://rdrr.io/r/base/print.html) for a compact overview,

- [`summary()`](https://rdrr.io/r/base/summary.html) for
  marginal/correlation diagnostics,

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for quick
  visualization of endpoint relationships.

## References

Cario, M. C., & Nelson, B. L. (1997). *Modeling and generating random
vectors with arbitrary marginal distributions and correlation matrix*
(pp. 1-19). Technical Report, Department of Industrial Engineering and
Management Sciences, Northwestern University, Evanston, Illinois.

## See also

[`endpoint_details`](https://boehringer-ingelheim.github.io/endpoints/reference/endpoint_details.md)
for endpoint specification details.

[`enrollment_details`](https://boehringer-ingelheim.github.io/endpoints/reference/enrollment_details.md)
for administrative censoring and stochastic enrollment options.

[`calibration_control`](https://boehringer-ingelheim.github.io/endpoints/reference/calibration_control.md)
for calibration tuning parameters.

[`corr_make`](https://boehringer-ingelheim.github.io/endpoints/reference/corr_make.md)
for creating correlation matrices.

See the
[`vignette("user_guide", package = "endpoints")`](https://boehringer-ingelheim.github.io/endpoints/articles/user_guide.md)
vignette for introductory examples, and the
[`vignette("longitudinal_data", package = "endpoints")`](https://boehringer-ingelheim.github.io/endpoints/articles/longitudinal_data.md)
vignette for simulating longitudinal data.

## Examples

``` r
library(endpoints)
## One continuous endpoint, two-arm trial
ep1 <- list(
  endpoint_type = "continuous",
  baseline_mean = 10,
  sd            = 2,
  trt_effect    = -1
)

sim1 <- makeData(
  correlation_matrix    = NULL,
  sample_size_per_group = 200,
  SEED                  = 1,
  endpoint_details      = list(ep1)
)

sim1
#> <makeDataSim>
#>   n =400
#>   n_arms =2
#>   endpoints =1
#>   target_correlation =FALSE
#>   single_endpoint_mode =TRUE
#>   control_only =FALSE
#>   n_by_arm = 0:200, 1:200
#>   endpoint_types = continuous
#> 
#>   data (head):
#>             Cont_1 trt
#> 1 8.74709238574752   0
#> 2  9.3475332695636   0
#> 3 10.3672866483698   0
#> 4 12.6595985189108   0
#> 5 8.32874275739157   0
#> ⋮                ⋮   ⋮
#> 
summary(sim1)
#> <summary.makeDataSim>
#> 
#> n_arms: 2 
#> 
#> Target correlation:
#> NULL
#> 
#> Estimated correlation (by arm):
#> 
#> arm_0:
#>        Cont_1
#> Cont_1      1
#> 
#> arm_1:
#>        Cont_1
#> Cont_1      1
#> 
#> Continuous endpoints (endpoint x arm):
#>   endpoint arm input_baseline_mean input_sd input_trt_effect est_baseline_mean
#> 1   Cont_1   0                  10        2                0          10.12674
#> 2   Cont_1   1                  10        2               -1          10.12674
#>   est_trt_effect est_resid_sd
#> 1        0.00000     1.773323
#> 2       -1.22838     1.920200

## Three correlated endpoints
R3 <- corr_make(
  num_endpoints = 3,
  values = rbind(
    c(1, 2, 0.20),
    c(1, 3, 0.10),
    c(2, 3, 0.15)
  )
)

ep_cont <- list(
  endpoint_type = "continuous",
  baseline_mean = 10,
  sd            = 2,
  trt_effect    = -1
)

ep_bin <- list(
  endpoint_type = "binary",
  baseline_prob = 0.30,
  trt_prob      = 0.45
)

ep_cnt <- list(
  endpoint_type = "count",
  baseline_mean = 8,
  trt_count     = 10,
  size          = 20,
  p_zero        = 0
)

sim3 <- makeData(
  correlation_matrix    = R3,
  sample_size_per_group = 1000,
  SEED                  = 123,
  endpoint_details      = list(ep_cont, ep_bin, ep_cnt),
  target_correlation    = TRUE
)

summary(sim3)
#> <summary.makeDataSim>
#> 
#> n_arms: 2 
#> 
#> Target correlation:
#>      [,1] [,2] [,3]
#> [1,]  1.0 0.20 0.10
#> [2,]  0.2 1.00 0.15
#> [3,]  0.1 0.15 1.00
#> 
#> Estimated correlation (by arm):
#> 
#> arm_0:
#>        Cont_1 Bin_1 Int_1
#> Cont_1  1.000 0.185 0.098
#> Bin_1   0.185 1.000 0.121
#> Int_1   0.098 0.121 1.000
#> 
#> arm_1:
#>        Cont_1 Bin_1 Int_1
#> Cont_1  1.000 0.164 0.046
#> Bin_1   0.164 1.000 0.171
#> Int_1   0.046 0.171 1.000
#> 
#> Continuous endpoints (endpoint x arm):
#>   endpoint arm input_baseline_mean input_sd input_trt_effect est_baseline_mean
#> 1   Cont_1   0                  10        2                0          9.906936
#> 2   Cont_1   1                  10        2               -1          9.906936
#>   est_trt_effect est_resid_sd
#> 1      0.0000000     2.033136
#> 2     -0.8925653     1.942604
#> 
#> Binary endpoints (endpoint x arm):
#>   endpoint arm input_baseline_prob input_trt_logOR input_trt_prob
#> 1    Bin_1   0                 0.3       0.0000000           0.30
#> 2    Bin_1   1                 0.3       0.6466272           0.45
#>   est_baseline_prob est_trt_logOR est_prob
#> 1              0.28     0.0000000    0.280
#> 2              0.28     0.7599398    0.454
#> 
#> Count endpoints (endpoint x arm):
#>   endpoint arm input_baseline_mean input_trt_logRR input_trt_mean input_size
#> 1    Int_1   0                   8       0.0000000              8         20
#> 2    Int_1   1                   8       0.2231436             10         20
#>   input_p_zero est_baseline_mean est_trt_logRR est_size obs_mean obs_p0
#> 1            0             8.037     0.0000000 20.49916    8.037  0.001
#> 2            0             8.037     0.2107994 20.49916    9.923  0.000

## One TTE endpoint with administrative censoring and exponential enrollment
ep_tte <- list(
  endpoint_type  = "tte",
  baseline_rate  = 1 / 24,
  trt_effect     = log(0.8),
  fatal_event    = TRUE
)

sim_tte <- makeData(
  correlation_matrix    = NULL,
  sample_size_per_group = 500,
  endpoint_details      = list(ep_tte),
  enrollment_details    = list(
    administrative_censoring    = 24,
    enrollment_distribution     = "exponential",
    enrollment_exponential_rate = 1 / 4
  )
)

summary(sim_tte)
#> <summary.makeDataSim>
#> 
#> n_arms: 2 
#> 
#> Target correlation:
#> NULL
#> 
#> Estimated correlation (by arm):
#> 
#> arm_0:
#>       TTE_1
#> TTE_1     1
#> 
#> arm_1:
#>       TTE_1
#> TTE_1     1
#> 
#> TTE endpoints (endpoint x arm):
#>   endpoint arm censor_col input_baseline_rate input_trt_logHR input_trt_HR
#> 1    TTE_1   0   Status_1          0.04166667       0.0000000          1.0
#> 2    TTE_1   1   Status_1          0.04166667      -0.2231436          0.8
#>   est_trt_logHR est_trt_HR obs_event_rate   exp_rate
#> 1     0.0000000  1.0000000          0.584 0.04455648
#> 2    -0.3631404  0.6954888          0.452 0.03085210
```
