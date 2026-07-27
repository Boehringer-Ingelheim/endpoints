# Specification format for `enrollment_details` used by `makeData()`

`enrollment_details` is a named list controlling stochastic subject
enrollment in
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md).

## Details

Enrollment times are represented on the trial calendar scale. That is,
`enrollTime` is the calendar time at which a subject enters the trial,
measured relative to the first randomized subject.

If omitted,
[`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
fills in defaults internally:


    list(
      enrollment_distribution        = "none",
      enrollment_exponential_rate    = NULL,
      piecewise_enrollment_cutpoints = NULL,
      piecewise_enrollment_rates     = NULL
    )

## Purpose

Enrollment determines when subjects enter the trial calendar. It is
distinct from subject-level follow-up limits and trial-level stopping
rules.

The trial calendar is generally constructed in three conceptual steps:

1.  enrollment is simulated using `enrollment_details`;

2.  the trial end time is determined using `trial_end_details`;

3.  subject-level available follow-up is derived using
    `followup_details`.

For a subject \\i\\, `enrollTime_i` is measured on the trial calendar
scale. These enrollment times can be used with any endpoint type. For
time-to-event endpoints, variables such as `TTE_1` are measured from
that subject's own enrollment or randomization time.

If the trial ends at calendar time \\\mathcal{T}\\, the subject's
available administrative follow-up is generally: \$\$ \max(0,
\mathcal{T} - T\_{E,i}), \$\$ where \\T\_{E,i}\\ is the subject's
enrollment time. If `followup_details$max_followup` is supplied,
available follow-up is further capped at that value.

If no enrollment model is supplied, `enrollment_distribution = "none"`,
all subjects are treated as enrolling at trial calendar time 0.

## Arguments

- `enrollment_distribution`:

  Character string giving the enrollment process. Must be one of:

  `"none"`

  :   All subjects enroll at time 0.

  `"exponential"`

  :   Homogeneous Poisson-process accrual. Inter-enrollment gaps are
      exponentially distributed with rate `enrollment_exponential_rate`.
      Requires `enrollment_exponential_rate`.

  `"piecewise"`

  :   Piecewise homogeneous Poisson-process accrual. Accrual rates vary
      by calendar interval. Requires `piecewise_enrollment_cutpoints`
      and `piecewise_enrollment_rates`.

- `enrollment_exponential_rate`:

  Numeric scalar greater than 0. Used only when
  `enrollment_distribution = "exponential"`.  
    
  This is the expected accrual rate per unit calendar time. For example,
  `enrollment_exponential_rate = 10` corresponds to an average of 10
  subjects enrolled per time unit. Inter-enrollment gaps are generated
  from an exponential distribution with this rate, and enrollment times
  are formed as cumulative sums of those gaps.

- `piecewise_enrollment_cutpoints`:

  Numeric vector of strictly increasing cutpoints defining calendar
  intervals for piecewise accrual. The first cutpoint should be 0.  
    
  For example, `c(0, 8, 16, 24)` defines three intervals: \\\[0,8)\\,
  \\\[8,16)\\, and \\\[16,24\]\\.

- `piecewise_enrollment_rates`:

  Numeric vector of non-negative accrual rates, one for each interval
  defined by `piecewise_enrollment_cutpoints`. Its length must be
  `length(piecewise_enrollment_cutpoints) - 1`.  
    
  At least one rate must be greater than 0. A rate of 0 can be used to
  encode an interval with no accrual. If the finite piecewise accrual
  window does not generate enough subjects,
  [`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
  continues enrollment beyond the final cutpoint using the final
  interval's rate. Therefore, if the final rate is 0 and the target
  sample size has not been reached, simulation will fail.

## No stochastic enrollment

When `enrollment_distribution = "none"`, all subjects are assigned:


    enrollTime = 0

In this case, if a trial end rule is used, all subjects share the same
calendar entry time. For example, if
`trial_end_details$type = "fixed_calendar"` and `trial_end_time = 80`,
each subject's calendar-derived available follow-up is 80, subject to
any `max_followup` cap.

## Exponential enrollment

When `enrollment_distribution = "exponential"`, enrollment follows a
homogeneous Poisson-process style accrual model.

The inter-enrollment gaps are generated as: \$\$ G_i \sim
\mathrm{Exponential}(\lambda), \$\$ where \\\lambda\\ is
`enrollment_exponential_rate`. Enrollment times are then generated as
cumulative sums: \$\$ T\_{E,i} = \sum\_{m=1}^{i} G_m. \$\$

The enrollment clock is shifted so the first enrolled subject has
`enrollTime = 0`. In
[`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md),
these calendar enrollment times are randomly assigned across generated
subject rows to avoid confounding treatment arm with enrollment order.

This model is appropriate when the expected enrollment rate is
approximately constant over calendar time.

## Piecewise enrollment

When `enrollment_distribution = "piecewise"`, enrollment follows a
piecewise homogeneous Poisson-process style accrual model.

Enrollment is simulated using exponential inter-arrival gaps. At any
point on the trial calendar, the waiting time to the next enrolled
subject is drawn using the accrual rate for the current calendar
interval: \$\$ G \sim \mathrm{Exponential}(\lambda_k), \$\$ where
\\\lambda_k\\ is the accrual rate for interval \\k\\.

If the proposed next enrollment time falls within the current interval,
that subject is enrolled and the process continues from the new
enrollment time. If the proposed next enrollment time crosses the next
cutpoint, no subject is enrolled before that boundary; the trial
calendar advances to the cutpoint, and the process continues using the
accrual rate for the next interval.

This construction supports start-up, ramp-up, plateau, slowdown, and
temporary enrollment pause patterns. Intervals with rate 0 are allowed
and represent periods with no accrual.

For example:


    enrollment_details <- list(
      enrollment_distribution        = "piecewise",
      piecewise_enrollment_cutpoints = c(0, 8, 24, 52),
      piecewise_enrollment_rates     = c(2, 8, 12)
    )

This represents expected accrual of:

- 2 subjects per time unit from 0 to 8,

- 8 subjects per time unit from 8 to 24,

- 12 subjects per time unit from 24 to 52.

If the requested sample size is not reached by the final cutpoint,
enrollment continues beyond the last cutpoint using the final interval's
rate. For example, with cutpoints `c(0, 12, 24, 52)` and rates
`c(1, 5, 10)`, the rate of 10 applies from 24 to 52 and also after 52 if
additional subjects are needed.

## Interaction with follow-up and trial-end rules

`enrollment_details` only determines when subjects enter the trial. It
can be used with continuous, binary, count, or TTE endpoints. By itself,
it does not alter endpoint values, and it does not by itself censor
time-to-event endpoints.

Administrative follow-up is determined by combining enrollment times
with `trial_end_details` and `followup_details`. For example:


    followup_details <- list(
      min_followup = 52,
      max_followup = 156
    )

    trial_end_details <- list(
      type = "last_patient_min_followup"
    )

Under this design, the trial ends when the last randomized subject has
at least 52 time units of follow-up: \$\$ \mathcal{T} =
\max_i(T\_{E,i}) + 52. \$\$

Subject-level available follow-up is then: \$\$ \min\\M, \max(0,
\mathcal{T} - T\_{E,i})\\, \$\$ where \\M\\ denotes `max_followup`.

For non-TTE endpoints, this same calendar logic still determines the
returned `enrollTime`, `availableFollowup`, and trial-calendar metadata,
even though no event time is administratively censored.

For a fixed-calendar trial:


    trial_end_details <- list(
      type = "fixed_calendar",
      trial_end_time = 80
    )

the trial ends at calendar time 80, and subject-level available
follow-up is: \$\$ \min\\M, \max(0, 80 - T\_{E,i})\\, \$\$ again with
\\M\\ denoting `max_followup`.

Subjects with `enrollTime` after the trial end time have
`availableFollowup = 0`. These rows may represent planned subjects who
would not contribute positive follow-up under the realized trial
calendar.

## Event-driven trials

For event-driven trials, enrollment times are used to convert
subject-level event times to calendar event times: \$\$
\mathrm{calendar\\ event\\ time}\_i = T\_{E,i} + T_i, \$\$ where \\T_i\\
is the subject-level time-to-event value.

For example:


    trial_end_details <- list(
      type = "event_driven",
      event_endpoint = "TTE_1",
      target_events = 150
    )

The trial end time is the calendar time at which the target number of
observed events is reached, subject to any additional options in
`trial_end_details`, such as `require_min_followup` or
`max_trial_duration`.

In event-driven trials,
[`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
returns the planned simulated sample. If the event target is reached
before all planned subjects have enrolled, some rows may have enrollment
times after the realized trial end time. These rows represent subjects
who would not contribute positive follow-up under the realized trial
calendar.

Users can define the analysis population using the realized trial end
time stored in the object metadata:


    dat <- as.data.frame(sim)
    trial_end <- sim$meta$trial_calendar$trial_end_time

    dat_analysis <- subset(dat, enrollTime < trial_end & availableFollowup > 0)

## Deprecated behavior

Earlier versions of
[`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
allowed `enrollment_details$administrative_censoring` and
`enrollment_distribution = "uniform"`.

These are no longer supported as part of the trial-calendar interface.
Administrative censoring is now represented through the combination of
`followup_details` and `trial_end_details`. Uniform enrollment has been
removed in favor of Poisson-process based accrual models.

If `administrative_censoring` is supplied in `enrollment_details`,
[`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
now errors and asks users to express administrative follow-up through
`followup_details` and `trial_end_details`.

## See also

[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
for the main simulation function.

[`endpoint_details`](https://boehringer-ingelheim.github.io/endpoints/reference/endpoint_details.md)
for endpoint specifications.

[`calibration_control`](https://boehringer-ingelheim.github.io/endpoints/reference/calibration_control.md)
for correlation-calibration settings.

## Examples

``` r
## No stochastic enrollment
enrollment_details <- list(
  enrollment_distribution = "none"
)

## Homogeneous Poisson-process enrollment
## Average accrual rate of 10 subjects per time unit
enrollment_details <- list(
  enrollment_distribution     = "exponential",
  enrollment_exponential_rate = 10
)

## Piecewise Poisson-process enrollment
enrollment_details <- list(
  enrollment_distribution        = "piecewise",
  piecewise_enrollment_cutpoints = c(0, 8, 24, 52),
  piecewise_enrollment_rates     = c(2, 8, 12)
)

## Example with last-patient-minimum-follow-up trial end
followup_details <- list(
  min_followup = 52,
  max_followup = 156
)

trial_end_details <- list(
  type = "last_patient_min_followup"
)

## Example with fixed calendar trial end
followup_details <- list(
  max_followup = 104
)

trial_end_details <- list(
  type = "fixed_calendar",
  trial_end_time = 80
)

## Example with event-driven trial end
trial_end_details <- list(
  type = "event_driven",
  event_endpoint = "TTE_1",
  target_events = 150
)
```
