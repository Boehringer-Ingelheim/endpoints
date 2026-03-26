# Specification format for `enrollment_details` used by `makeData()`

`enrollment_details` is a named list controlling administrative
censoring and stochastic enrollment in
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md).

If omitted,
[`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
fills in defaults internally:

    list(
      administrative_censoring       = NULL,
      enrollment_distribution        = "none",
      enrollment_exponential_rate    = NULL,
      piecewise_enrollment_cutpoints = NULL,
      piecewise_enrollment_rates     = NULL
    )

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

- `administrative_censoring`:

  Numeric scalar or `NULL`. If non-`NULL` and positive, this defines the
  administrative end of follow-up. Any TTE outcome occurring after this
  limit is administratively censored.

- `enrollment_distribution`:

  Character string giving the enrollment-time distribution. Must be one
  of:

  `"none"`

  :   All subjects enroll at time 0.

  `"uniform"`

  :   Enrollment times are drawn from \\U(0,\mathcal{A})\\.

  `"exponential"`

  :   Enrollment times are drawn from an exponential distribution,
      truncated at \\\mathcal{A}\\. Requires
      `enrollment_exponential_rate`.

  `"piecewise"`

  :   Enrollment times are generated from a piecewise exponential model
      over user-specified intervals. Requires
      `piecewise_enrollment_cutpoints` and `piecewise_enrollment_rates`.

- `enrollment_exponential_rate`:

  Numeric scalar (\> 0). Used only when
  `enrollment_distribution = "exponential"`. Gives the exponential rate
  governing enrollment times.

- `piecewise_enrollment_cutpoints`:

  Numeric vector of strictly increasing cutpoints defining the intervals
  for piecewise enrollment. For example, `c(0, 8, 16, 24)` defines three
  intervals: \[0,8), \[8,16), and \[16,24\].

- `piecewise_enrollment_rates`:

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

## Examples

``` r
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
