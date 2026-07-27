# Specification format for `trial_end_details` used by `makeData()`

`trial_end_details` is a named list controlling when the trial ends on
the calendar scale in
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md).

## Details

If omitted,
[`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
fills in defaults internally:


    list(
      type                 = "none",
      trial_end_time       = NULL,
      event_endpoint       = NULL,
      target_events        = NULL,
      require_min_followup = FALSE,
      max_trial_duration   = NULL,
      target_not_reached   = "error"
    )

## Purpose

`trial_end_details` determines the realized trial end time on the
calendar scale. This is the time at which administrative follow-up stops
and final analysis is assumed to occur.

It works together with:

- `enrollment_details`, which determines when each subject enters the
  trial, and

- `followup_details`, which determines subject-level follow-up limits
  once the trial end time is known.

These trial-calendar features can be used with continuous, binary,
count, or time-to-event endpoints. For non-TTE endpoints, the main
effect is on the returned calendar quantities such as `enrollTime`,
`availableFollowup`, and `sim$meta$trial_calendar`.

For time-to-event endpoints, enrollment times and event times also
combine to form calendar event times: \$\$ \mathrm{calendar\\ event\\
time}\_i = T\_{E,i} + T_i, \$\$ where \\T\_{E,i}\\ is enrollment time
and \\T_i\\ is subject-level follow-up time to event or censoring.

## Arguments

- `type`:

  Character string specifying the trial-ending rule. Must be one of:

  `"none"`

  :   No trial-level stopping rule is imposed. The realized trial end
      time is treated as infinite, subject only to any subject-level
      `max_followup` cap.

  `"fixed_calendar"`

  :   End the trial at a user-specified calendar time given by
      `trial_end_time`.

  `"last_patient_min_followup"`

  :   End the trial when the last randomized subject has reached
      `followup_details$min_followup`.

  `"event_driven"`

  :   End the trial when a target number of observed events has accrued
      on a chosen TTE endpoint. This option is specific to time-to-event
      endpoints and may be combined with further restrictions such as
      `require_min_followup` or `max_trial_duration`.

- `trial_end_time`:

  Non-negative numeric scalar used when `type = "fixed_calendar"`.

- `event_endpoint`:

  Endpoint used to drive event-driven stopping. This may be supplied as
  a TTE ordinal such as `1`, or as a character label such as `"TTE_1"`,
  `"Status_1"`, or the raw endpoint index format `"V3"`, provided it
  resolves to a TTE endpoint.

- `target_events`:

  Positive integer giving the number of observed events required for
  `type = "event_driven"`.

- `require_min_followup`:

  Logical scalar used in event-driven trials. If `TRUE`, the trial
  cannot end before the last randomized subject has reached
  `followup_details$min_followup`, even if the event target is reached
  earlier.

- `max_trial_duration`:

  Optional positive numeric scalar giving a hard upper bound on calendar
  trial duration. This can be used either as a global cap or as a
  fallback when an event target is not reached.

- `target_not_reached`:

  Character string controlling what happens when an event-driven design
  does not reach `target_events`. Must be one of:

  `"error"`

  :   Stop with an error.

  `"max_trial_duration"`

  :   End the trial at `max_trial_duration`. Requires
      `max_trial_duration`.

  `"last_patient_min_followup"`

  :   End the trial when the last randomized subject reaches
      `followup_details$min_followup`. Requires
      `followup_details$min_followup`.

## Fixed calendar trial end

For a fixed-calendar design:


    trial_end_details <- list(
      type = "fixed_calendar",
      trial_end_time = 36
    )

the trial ends at calendar time 36, regardless of how many subjects have
enrolled or how many events have occurred by then.

This design is often useful when the database lock or final analysis
date is prespecified on the calendar scale, including for non-TTE
endpoint simulations where users still want realistic enrollment timing
and available-follow-up metadata.

## Last-patient-minimum-follow-up trial end

For a last-patient-minimum-follow-up design:


    followup_details <- list(
      min_followup = 24
    )

    trial_end_details <- list(
      type = "last_patient_min_followup"
    )

the realized trial end time is: \$\$ \mathcal{T} = \max_i(T\_{E,i}) +
24. \$\$

This design is useful when the analysis should not occur until the
latest enrolled subject has had enough time under observation. It can be
used for any endpoint type when users want the trial calendar to reflect
a minimum amount of follow-up for late enrollees.

## Event-driven trial end

For an event-driven design:


    trial_end_details <- list(
      type = "event_driven",
      event_endpoint = "TTE_1",
      target_events = 150
    )

the trial end time is the calendar time at which the chosen endpoint
first accumulates the requested number of observed events, subject to
any `max_followup`, `require_min_followup`, and `max_trial_duration`
restrictions.

This option is only available when at least one TTE endpoint is present.
It is especially relevant when enrollment is staggered, because the
calendar time of an event depends on both `enrollTime` and the
subject-level event time.

## If the event target is not reached

In event-driven designs, the requested target may be unattainable under
the simulated event rate, censoring rate, follow-up window, or sample
size. The `target_not_reached` argument controls how to handle this
case.

For example:


    trial_end_details <- list(
      type = "event_driven",
      event_endpoint = "TTE_1",
      target_events = 500,
      target_not_reached = "max_trial_duration",
      max_trial_duration = 48
    )

or


    followup_details <- list(
      min_followup = 18
    )

    trial_end_details <- list(
      type = "event_driven",
      event_endpoint = "TTE_1",
      target_events = 500,
      target_not_reached = "last_patient_min_followup"
    )

## Trial-calendar metadata

[`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
stores realized calendar outputs in:


    sim$meta$trial_calendar

This metadata includes the realized trial end time, the stopping reason,
whether an event target was reached, and the number of events observed
at the realized trial end.

For non-TTE endpoint simulations, this metadata still records the
realized calendar design even though no event-driven event count is
involved.

In event-driven designs,
[`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
returns the planned sample. If the realized event-driven stop occurs
before all planned subjects have accrued positive follow-up, some rows
may have `availableFollowup = 0`.

## See also

[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
for the main simulation function.

[`followup_details`](https://boehringer-ingelheim.github.io/endpoints/reference/followup_details.md)
for subject-level follow-up rules.

[`enrollment_details`](https://boehringer-ingelheim.github.io/endpoints/reference/enrollment_details.md)
for enrollment-time settings.

## Examples

``` r
## Fixed calendar analysis
trial_end_details <- list(
  type = "fixed_calendar",
  trial_end_time = 36
)

## Last-patient-minimum-follow-up design
followup_details <- list(
  min_followup = 24
)

trial_end_details <- list(
  type = "last_patient_min_followup"
)

## Event-driven design
trial_end_details <- list(
  type = "event_driven",
  event_endpoint = "TTE_1",
  target_events = 50
)

## Event-driven fallback to a maximum trial duration
trial_end_details <- list(
  type = "event_driven",
  event_endpoint = "TTE_1",
  target_events = 500,
  target_not_reached = "max_trial_duration",
  max_trial_duration = 48
)

## Small makeData() example
ep_tte <- list(
  endpoint_type = "tte",
  baseline_rate = rate_from_prob(
    target_prob = 0.25,
    mode = "admin",
    admin_time = 12
  ),
  trt_effect = log(0.80),
  censoring_rate = 0.01
)

sim <- makeData(
  correlation_matrix = NULL,
  sample_size_per_group = 12,
  SEED = 2,
  endpoint_details = list(ep_tte),
  enrollment_details = list(
    enrollment_distribution = "exponential",
    enrollment_exponential_rate = 5
  ),
  followup_details = list(
    max_followup = 24
  ),
  trial_end_details = list(
    type = "event_driven",
    event_endpoint = "TTE_1",
    target_events = 6
  )
)

sim$meta$trial_calendar
#> $trial_end_time
#> [1] 13.89784
#> 
#> $trial_end_reason
#> [1] "event_target"
#> 
#> $event_target_reached
#> [1] TRUE
#> 
#> $n_events_at_end
#> [1] 6
#> 
```
