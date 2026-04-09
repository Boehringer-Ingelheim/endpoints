# Specification format for `endpoint_details` used by `makeData()`

`endpoint_details` is the main endpoint-specification input to
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md).
It is a non-empty list in which each element is itself a named list
describing one endpoint.

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

## General structure

For \\p\\ endpoints, `endpoint_details` should be a list of length
\\p\\:

    endpoint_details = list(
      ep1, ep2, ..., ep_p
    )

Each `ep_j` is a named list containing the parameters for endpoint
\\j\\.

## Continuous endpoints

Continuous endpoints are generated using the Gaussian distribution.

Fields:

- endpoint_type:

  Set to `"continuous"` or a recognized alias such as `"normal"`.

- baseline_mean:

  Numeric scalar giving the control-group mean.

- sd:

  Numeric scalar or numeric vector. If scalar, a common SD is used
  across all arms. If vector, it must have length equal to the total
  number of arms.

- trt_effect:

  Optional. Numeric scalar or numeric vector giving treatment-group mean
  shifts relative to control. A scalar is used for 2-arm trials;
  otherwise the length must be `K - 1`, where `K` is the total number of
  arms. If `NULL`, only data for the control group is generated. If
  `NULL` and `arm_mode = "full"`, `trt_effect` is taken to be 0 for all
  active arms.

## Binary endpoints

For binary outcomes, the event probability is modeled on the logit
scale.

Fields:

- endpoint_type:

  Set to `"binary"`.

- baseline_prob:

  Numeric scalar in `(0, 1)` giving the control-group event probability.

- trt_prob:

  Optional. Numeric scalar or vector giving treatment-group event
  probabilities directly.

- trt_effect:

  Optional. Numeric scalar or vector giving treatment-group effects on
  the log-odds scale.

Exactly one of `trt_prob` or `trt_effect` may be supplied. If neither is
supplied, a control-only binary endpoint is generated, or a no-effect
endpoint in a multi-arm setting depending on `arm_mode`.

## Count endpoints

Count endpoints are generated using the negative binomial or
zero-inflated negative binomial distribution.

Required fields:

- endpoint_type:

  Set to `"count"`.

- baseline_mean:

  Numeric scalar (\> 0) giving the control-group mean count.

- size:

  Numeric scalar (\> 0) controlling the negative-binomial dispersion.
  Larger values reduce overdispersion. Set to a large value to
  approximate the Poisson distribution.

- trt_count:

  Optional. Numeric scalar or vector giving treatment-group mean counts
  directly.

- trt_effect:

  Optional. Numeric scalar or vector giving treatment-group effects on
  the log rate-ratio scale.

- p_zero:

  Optional. Numeric scalar in \\\[0,1\]\\ giving the structural-zero
  probability. If omitted, defaults to 0.

Exactly one of `trt_count` or `trt_effect` may be supplied. In the
parameterization used, the variance is \\\mu + \frac{\mu^2}{\phi}\\,
where `size` controls \\\phi\\. The `p_zero` parameter is independent of
the copula.

## Time-to-event endpoints

Time-to-event endpoints are generated from exponential event-time
distributions.

Required fields:

- endpoint_type:

  Set to `"time-to-event"` or a recognized alias such as `"tte"`.

- baseline_rate:

  Numeric scalar (\> 0) giving the control-group exponential event-rate
  parameter \\\lambda\\. The mean event time is therefore \\1/\lambda\\.

Optional fields:

- trt_effect:

  Numeric scalar or vector giving treatment-group effects on the log
  hazard-ratio scale.

- censoring_rate:

  Numeric scalar (\> 0) giving the independent exponential
  censoring-rate parameter. If `NULL`, no random censoring is applied.

- fatal_event:

  Logical scalar. If `TRUE`, the endpoint is treated as fatal and
  censors later TTE endpoints according to the package's fatal/non-fatal
  rules.

## Arm-specific lengths

For active-arm parameters:

- scalar values should be used for two-arm trials,

- vectors must generally have length `K - 1`, where `K` is the total
  number of arms,

- for continuous `sd`, vectors must have length `K`, since the
  control-group SD is also included explicitly.

The number of arms is determined from the supplied endpoint
specifications. Optionally, users may use `arm_mode` in
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md),
but this option should be used with care.

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

  `Status_1`, `Status_2`, ..., where `1 = event` and `0 = censored`

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
```
