# Specification format for `followup_details` used by `makeData()`

`followup_details` is a named list controlling subject-level
administrative follow-up in
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md).

## Details

If omitted,
[`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
fills in defaults internally:


    list(
      min_followup = NULL,
      max_followup = NULL
    )

## Purpose

`followup_details` determines how much calendar follow-up each subject
can contribute once enrollment and the trial end time have been
determined.

It is distinct from:

- `enrollment_details`, which determines when subjects enter the trial
  calendar, and

- `trial_end_details`, which determines when the trial stops on the
  calendar scale.

These three components define the realized trial calendar for any
endpoint type:

1.  `enrollment_details` generates `enrollTime`;

2.  `trial_end_details` determines the realized trial end time
    \\\mathcal{T}\\;

3.  `followup_details` determines each subject's `availableFollowup`.

For time-to-event endpoints, this same calendar structure also
determines administrative censoring of event times and statuses.

A subject enrolled at calendar time \\T\_{E,i}\\ can contribute at most:
\$\$ \max(0, \mathcal{T} - T\_{E,i}) \$\$ units of calendar follow-up
before the trial ends. If `max_followup` is supplied, this is further
capped at: \$\$ \min\\M, \max(0, \mathcal{T} - T\_{E,i})\\, \$\$ where
\\M\\ denotes `max_followup`.

## Arguments

- `min_followup`:

  Optional non-negative numeric scalar. This is the minimum follow-up
  duration that certain trial-end rules may require for the last
  randomized subject.  
    
  By itself, `min_followup` does not censor subjects and does not force
  the trial to continue. It becomes operational when used together with
  a trial-end rule such as
  `trial_end_details$type = "last_patient_min_followup"` or
  `trial_end_details$require_min_followup = TRUE` in an event-driven
  trial.

- `max_followup`:

  Optional positive numeric scalar giving the maximum observable
  follow-up for any subject.  
    
  This affects the realized trial calendar for any endpoint type by
  capping `availableFollowup`. For time-to-event endpoints, it also acts
  as a subject-level administrative censoring cap. If a subject would
  otherwise remain under observation longer than `max_followup`, their
  observed TTE time is truncated at `max_followup` and the corresponding
  event indicator is administratively censored when appropriate. For
  non-TTE endpoints, `max_followup` changes the returned calendar
  metadata but does not modify the endpoint values themselves.

## Common use cases

- Only `max_followup`:

  Use this when each subject can be followed for at most a fixed amount
  of time, regardless of how long the trial remains open.

- Only `min_followup`:

  Use this when the design requires the last randomized subject to have
  at least a certain amount of follow-up before final analysis, for
  example with `trial_end_details$type = "last_patient_min_followup"`.

- Both `min_followup` and `max_followup`:

  Use this when the design requires a minimum observation window for
  late enrollees, while also capping how long early enrollees can remain
  under observation.

## Interaction with trial-end rules

The effect of `followup_details` depends on the trial-ending rule.

For a fixed-calendar trial:


    trial_end_details <- list(
      type = "fixed_calendar",
      trial_end_time = 36
    )

subjects enrolled earlier will usually have more potential follow-up
than subjects enrolled later, but each subject is still capped at
`max_followup` if it is supplied.

For a last-patient-minimum-follow-up design:


    followup_details <- list(
      min_followup = 24,
      max_followup = 36
    )

    trial_end_details <- list(
      type = "last_patient_min_followup"
    )

the realized trial end is: \$\$ \mathcal{T} = \max_i(T\_{E,i}) + 24.
\$\$

For non-TTE endpoints, these rules still determine `enrollTime`,
`availableFollowup`, and the metadata in `sim$meta$trial_calendar`, even
though the endpoint values themselves are not administratively censored.

For an event-driven design, `max_followup` determines whether late TTE
events remain observable, and `min_followup` can optionally be enforced
through `trial_end_details$require_min_followup = TRUE`. Event-driven
trial end is specific to TTE endpoints.

## Returned data

When a trial-calendar feature is active,
[`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
returns an `availableFollowup` column in the simulated dataset. This
stores the realized administrative follow-up available to each subject
under the combination of enrollment timing, trial end time, and any
`max_followup` cap.

Subjects with `availableFollowup == 0` are subjects whose planned
enrollment occurs at or after the realized trial end time. They remain
in the simulated planned dataset but typically do not contribute to the
realized analysis population.

## See also

[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
for the main simulation function.

[`enrollment_details`](https://boehringer-ingelheim.github.io/endpoints/reference/enrollment_details.md)
for enrollment-time settings.

[`trial_end_details`](https://boehringer-ingelheim.github.io/endpoints/reference/trial_end_details.md)
for calendar stopping rules.

## Examples

``` r
## Follow each subject for at most 36 time units
followup_details <- list(
  max_followup = 36
)

## Require at least 24 time units of follow-up for the
## last randomized subject
followup_details <- list(
  min_followup = 24
)

## Combine a minimum and maximum follow-up rule
followup_details <- list(
  min_followup = 24,
  max_followup = 36
)

## Small makeData() example
ep_tte <- list(
  endpoint_type = "tte",
  baseline_rate = rate_from_prob(
    target_prob = 0.30,
    mode = "admin",
    admin_time = 12
  ),
  trt_effect = log(0.70),
  censoring_rate = 0.01
)

sim <- makeData(
  correlation_matrix = NULL,
  sample_size_per_group = 10,
  SEED = 1,
  endpoint_details = list(ep_tte),
  enrollment_details = list(
    enrollment_distribution = "exponential",
    enrollment_exponential_rate = 4
  ),
  followup_details = list(
    min_followup = 12,
    max_followup = 18
  ),
  trial_end_details = list(
    type = "last_patient_min_followup"
  )
)

dat <- as.data.frame(sim)
head(dat[, c("TTE_1", "Status_1", "enrollTime", "availableFollowup")])
#>       TTE_1 Status_1 enrollTime availableFollowup
#> 1 10.381792        1   0.449141          15.21788
#> 2 13.550787        0   2.116230          13.55079
#> 3 14.505765        0   1.161252          14.50576
#> 4 12.775784        0   2.891233          12.77578
#> 5  7.578268        1   1.342556          14.32446
#> 6 14.581200        0   1.085817          14.58120
```
