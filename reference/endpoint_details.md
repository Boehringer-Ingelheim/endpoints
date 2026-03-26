# Specification format for `endpoint_details` used by `makeData()`

`endpoint_details` is the main endpoint-specification input to
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md).
It is a non-empty list in which each element is itself a named list
describing one endpoint.

`enrollment_details` is a named list controlling administrative
censoring and stochastic enrollment in
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md).

## Details

Each endpoint specification must include `endpoint_type`, and then
additional fields depending on the endpoint family. Supported endpoint
types are:

- continuous,

- binary,

- count,

- time-to-event.

Aliases such as `"normal"`, `"gaussian"`, `"bin"`, `"nb"`, `"zinb"`, and
`"tte"` may also be accepted by the package's internal normalization
routines.

If omitted,
[`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
fills in defaults internally:

    list( administrative_censoring = NULL, enrollment_distribution = "none",
    enrollment_exponential_rate = NULL, piecewise_enrollment_cutpoints = NULL,
    piecewise_enrollment_rates = NULL ) 

## General structure

For \\p\\ endpoints, `endpoint_details` should be a list of length
\\p\\:

     endpoint_details = list(
    ep1, ep2, ..., ep_p ) 

Each `ep_j` is a named list containing the parameters for endpoint
\\j\\.

## Continuous endpoints

Continuous endpoints are generated using the Gaussian distribution.

Fields:

- list("endpoint_type"):

  Set to `"continuous"` (or a recognized alias such as `"normal"`).

- list("baseline_mean"):

  Numeric scalar giving the control-group mean.

- list("sd"):

  Numeric scalar or numeric vector. If scalar, a common SD is used
  across all arms. If vector, it must have length equal to the total
  number of arms.

- list("trt_effect"):

  Optional. Numeric scalar or numeric vector giving treatment-group mean
  shifts relative to control. A scalar is used for 2-arm trials;
  otherwise the length must be `K - 1`, where `K` is the total number of
  arms. If `NULL`, only data for the control group is generated. If
  `NULL` and `arm_mode = "full"`, trt_effect = 0 for all arms.

## Binary endpoints

For binary outcomes, the event probability is modeled on the logit
scale.

Fields:

- list("endpoint_type"):

  Set to `"binary"`.

- list("baseline_prob"):

  Numeric scalar in `(0, 1)` giving the control-group event probability.

- list("trt_prob"):

  Optional. Numeric scalar or vector giving treatment-group event
  probabilities directly.

- list("trt_effect"):

  Optional. Numeric scalar or vector giving treatment-group effects on
  the log-odds scale.

Exactly one of `trt_prob` or `trt_effect` may be supplied. If neither is
supplied, a control-only binary endpoint is generated (or a no-effect
endpoint in a multi-arm setting, depending on `arm_mode`).

## Count endpoints

Count endpoints are generated using the negative binomial or
zero-inflated negative binomial distribution.

Required fields:

- list("endpoint_type"):

  Set to `"count"`.

- list("baseline_mean"):

  Numeric scalar (\> 0) giving the control-group mean count.

- list("size"):

  Numeric scalar (\> 0) controlling the negative-binomial dispersion.
  Larger values reduce overdispersion. Set to a large value to
  approximate the Poisson distribution.

- list("trt_count"):

  Optional. Numeric scalar or vector giving treatment-group mean counts
  directly.

- list("trt_effect"):

  Optional. Numeric scalar or vector giving treatment-group effects on
  the log rate-ratio scale.

- list("p_zero"):

  Optional. Numeric scalar in \\\[0,1\]\\ giving the structural-zero
  probability. If omitted, defaults to 0.

Exactly one of `trt_count` or `trt_effect` may be supplied. In the
parameterization used, the variance of the distribution is \\\mu +
\frac{\mu^2}{\phi}\\, `size` controls \\\phi\\. The `p_zero` parameter
is indpendent of the copula.

## Time-to-event endpoints

Time-to-event endpoints are generated from exponential event-time
distributions.

Required fields:

- list("endpoint_type"):

  Set to `"time-to-event"` or a recognized alias such as `"tte"`.

- list("baseline_rate"):

  Numeric scalar (\>0) giving the control-group exponential event-rate
  parameter \\\lambda\\. The mean event time is therefore \\1/\lambda\\.

Optional fields:

- list("trt_effect"):

  Numeric scalar or vector giving treatment-group effects on the log
  hazard-ratio scale.

- list("censoring_rate"):

  Numeric scalar (\>0) giving the independent exponential censoring-rate
  parameter. If `NULL`, no random censoring.

- list("fatal_event"):

  Logical scalar. If `TRUE`, the endpoint is treated as fatal and
  censors later TTE endpoints according to the package's fatal/non-fatal
  rules.

## Arm-specific lengths

For active-arm parameters:

- scalar values should be used for two-arm trials;

- vectors must generally have length `K - 1`, where `K` is the total
  number of arms;

- for continuous `sd`, vectors must have length `K`, since the
  control-group SD is also included explicitly.

The number of arms is determined from the supplied endpoint
specifications. Optionally, users may use `arm_mode` in
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)`, but this option should be used with care.`.

## Internal endpoint naming

In the returned dataset, endpoints are automatically renamed by type:

- Continuous:

  `Cont_1`, `Cont_2`, ...

- Binary:

  `Bin_1`, `Bin_2`, ...

- Count:

  `Int_1`, `Int_2`, ...

- Time-to-event:

  `TTE_1`, `TTE_2`, ...

- TTE event indicators:

  `Status_1`, `Status_2`, ... (`1 = event`, `0 = censored`)

## Purpose

Administrative censoring and stochastic enrollment affect *observed
follow-up*. This primarily affects time-to-event endpoints, but may also
be of use for possible interim analyses.

If no enrollment model is supplied (`enrollment_distribution = "none"`),
all subjects are treated as enrolling at time 0.

If enrollment is stochastic and administrative censoring is enabled,
each subject receives an `enrollTime`, and maximum available follow-up
is reduced from \\\mathcal{A}\\ to \\\mathcal{A} - T_E\\, where
\\\mathcal{A}\\ is the administrative censoring time and \\T_E\\ is the
subject's enrollment time.

## Arguments

- list("administrative_censoring"):

  Numeric scalar or `NULL`. If non-`NULL` and positive, this defines the
  administrative end of follow-up. Any TTE outcome occurring after this
  limit is administratively censored.

- list("enrollment_distribution"):

  Character string giving the enrollment-time distribution. Must be one
  of:

  list("\\none\\")

  :   All subjects enroll at time 0.

  list("\\uniform\\")

  :   Enrollment times are drawn from \\U(0,\mathcal{A})\\.

  list("\\exponential\\")

  :   Enrollment times are drawn from an exponential distribution,
      truncated at \\\mathcal{A}\\. Requires
      `enrollment_exponential_rate`.

  list("\\piecewise\\")

  :   Enrollment times are generated from a piecewise exponential model
      over user-specified intervals. Requires
      `piecewise_enrollment_cutpoints` and `piecewise_enrollment_rates`.

- list("enrollment_exponential_rate"):

  Numeric scalar (\> 0). Used only when
  `enrollment_distribution = "exponential"`. Gives the exponential rate
  governing enrollment times.

- list("piecewise_enrollment_cutpoints"):

  Numeric vector of strictly increasing cutpoints defining the intervals
  for piecewise enrollment. For example, `c(0, 8, 16, 24)` defines three
  intervals: `[0,8)`, `[8,16)`, and `[16,24]`.

- list("piecewise_enrollment_rates"):

  Numeric vector of positive exponential rates, one for each interval
  defined by `piecewise_enrollment_cutpoints`. Its length must be
  `length(piecewise_enrollment_cutpoints) - 1`.

## Administrative censoring only

If `administrative_censoring` is supplied and
`enrollment_distribution = "none"`, all subjects are assumed to enter at
time 0 and remain under observation until the administrative end of
follow-up (unless they experience an event or are independently censored
earlier).

## Uniform enrollment

Under uniform enrollment, subjects enter uniformly over the interval
\\\[0,\mathcal{A}\]\\. This induces heterogeneous follow-up, with later
enrollees having less potential follow-up time before administrative
censoring.

## Exponential enrollment

Under exponential enrollment, subject entry times are sampled from an
exponential distribution and truncated at the administrative censoring
time. This can be used to represent rapid early enrollment followed by a
tapering accrual pattern.

## Piecewise exponential enrollment

Under piecewise enrollment, the waiting time to enrollment is simulated
interval-by-interval. Within each interval, an exponential waiting time
is drawn using that interval's rate. If the draw exceeds the remaining
interval width, the process moves to the next interval.

This allows users to encode enrollment patterns such as slow-then-fast
accrual, fast-then-slow accrual, or more complex piecewise shapes. See
the
[`vignette("user_guide", package = "endpoints")`](https://boehringer-ingelheim.github.io/endpoints/articles/user_guide.md)
vignette for an example of how to use this option to achieve certain
enrollment proportions in certain enrollment periods.

## See also

[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
for the main simulation routine.

[`enrollment_details`](https://boehringer-ingelheim.github.io/endpoints/reference/enrollment_details.md)
for follow-up and enrollment settings.

[`calibration_control`](https://boehringer-ingelheim.github.io/endpoints/reference/calibration_control.md)
for correlation-calibration controls.

## Examples

``` r
## Continuous endpoint
ep_cont <- list(
  endpoint_type = "continuous",
  baseline_mean = 10,
  sd            = c(2, 3),
  trt_effect    = -1
)

## Binary endpoint
ep_bin <- list(
  endpoint_type = "binary",
  baseline_prob = 0.30,
  trt_prob      = 0.45
)

## Count endpoint
ep_cnt <- list(
  endpoint_type = "count",
  baseline_mean = 8,
  trt_count     = 10,
  size          = 20,
  p_zero        = 0.10
)

## Time-to-event endpoint
ep_tte <- list(
  endpoint_type  = "tte",
  baseline_rate  = 1 / 24,
  trt_effect     = log(0.80),
  censoring_rate = 1 / 216,
  fatal_event    = TRUE
)

endpoint_details <- list(ep_cont, ep_bin, ep_cnt, ep_tte)


## Administrative censoring only
enrollment_details <- list(
  administrative_censoring = 24
)

## Exponential enrollment
enrollment_details <- list(
  administrative_censoring    = 24,
  enrollment_distribution     = "exponential",
  enrollment_exponential_rate = 1 / 4
)
head(enrollment_details$data)
#> NULL

## Piecewise exponential enrollment
enrollment_details <- list(
  administrative_censoring       = 24,
  enrollment_distribution        = "piecewise",
  piecewise_enrollment_cutpoints = c(0, 8, 16, 24),
  piecewise_enrollment_rates     = c(0.02, 0.06, 0.50)
)
```
