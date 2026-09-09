# Changelog

## endpoints 0.1.3

- Time-to-event endpoints now support arm-specific independent random
  censoring rates. Supply `censoring_rate` as a length-`K` vector
  ordered by treatment arm; a scalar value continues to apply the same
  censoring rate to all arms.

## endpoints 0.1.2

- Added stochastic enrollment with homogeneous exponential and piecewise
  Poisson-process accrual through `enrollment_details`.
- Added subject-level follow-up limits through `followup_details` and
  calendar trial-ending rules through `trial_end_details`, including
  fixed-calendar, last-patient-minimum-follow-up, and event-driven
  designs.
- Piecewise enrollment now supports an open-ended final accrual interval
  by supplying `piecewise_enrollment_rates` with the same length as
  `piecewise_enrollment_cutpoints`.
- [`plot.makeDataSim()`](https://boehringer-ingelheim.github.io/endpoints/reference/plot.makeDataSim.md)
  no longer depends on `GGally` and now draws its plot matrix directly.
  As a result, it returns `invisible(NULL)` rather than a ggplot object,
  so plots can no longer be customized afterward with ggplot2’s `+`
  syntax.
- [`plot.makeDataSim()`](https://boehringer-ingelheim.github.io/endpoints/reference/plot.makeDataSim.md)
  gains a `title` argument for overriding the default `"Arm <k>"` plot
  title.
- `enrollment_details$administrative_censoring` is no longer supported.
  Use `followup_details` and `trial_end_details` to express
  administrative follow-up instead.
- `enrollment_distribution = "uniform"` is no longer supported. Use
  `"none"`, `"exponential"`, or `"piecewise"` enrollment instead.

## endpoints 0.1.1

- Added package website configuration and package website links.

## endpoints 0.1.0

- Initial release.
