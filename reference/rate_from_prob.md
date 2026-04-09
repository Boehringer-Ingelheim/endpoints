# Solve exponential rate parameters from target observed event probabilities

`rate_from_prob()` is a helper for converting a user-specified observed
event probability into an exponential rate parameter under several
common time-to-event design settings.

## Usage

``` r
rate_from_prob(
  target_prob,
  mode = c("simple", "admin", "semi-competing"),
  event_rate = NULL,
  admin_time = NULL,
  fatal_event_rate = NULL,
  fatal_censor_rate = NULL,
  nonfatal_event_rate = NULL
)
```

## Arguments

- target_prob:

  A numeric scalar in `(0, 1)` giving the target observed event
  probability.

  The meaning depends on `mode`:

  - for `"simple"`, the desired probability that the event is observed
    before independent random censoring;

  - for `"admin"`, the desired probability that the event occurs before
    the administrative follow-up limit;

  - for `"semi-competing"`, the desired approximate probability that the
    secondary non-fatal event is observed first.

- mode:

  Character string indicating which probability-to-rate relationship to
  use. Must be one of:

  simple

  :   One time-to-event endpoint with independent exponential censoring.

  admin

  :   One time-to-event endpoint with administrative censoring only.

  semi-competing

  :   Approximate semi-competing risks calculation for a secondary
      non-fatal endpoint.

- event_rate:

  Numeric scalar giving the exponential event rate \\\lambda_e\\.

  Used only when `mode = "simple"`.

- admin_time:

  Numeric scalar giving the administrative censoring time horizon \\A\\.

  Used only when `mode = "admin"`.

- fatal_event_rate:

  Numeric scalar giving the event rate for the primary fatal event,
  denoted \\\lambda\_{e1}\\.

  Used only when `mode = "semi-competing"`.

- fatal_censor_rate:

  Numeric scalar giving the independent censoring rate for the fatal
  endpoint, denoted \\\lambda\_{c1}\\.

  Used only when `mode = "semi-competing"`.

- nonfatal_event_rate:

  Numeric scalar giving the event rate for the secondary non-fatal
  endpoint, denoted \\\lambda\_{e2}\\.

  Used only when `mode = "semi-competing"`.

## Value

A numeric scalar giving the solved rate parameter.

Depending on `mode`, this is:

- the censoring rate (`"simple"`),

- the event rate (`"admin"`), or

- the approximate non-fatal censoring rate (`"semi-competing"`).

## Details

Depending on `mode`, the function solves for:

- an independent exponential censoring rate (`"simple"`),

- an exponential event rate under pure administrative censoring
  (`"admin"`), or

- an approximate censoring rate for a secondary non-fatal endpoint in a
  semi-competing risks setting (`"semi-competing"`).

This is intended as a design-stage helper for choosing rate parameters
before calling
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md).
Not all cases are covered. In more complicated settings, iteration may
be necessary.

## Mode = `"simple"`

Assume a single event time \\T \sim \mathrm{Exp}(\lambda_e)\\ and an
independent censoring time \\C \sim \mathrm{Exp}(\lambda_c)\\.

The probability of observing the event is

\$\$ P(T \< C) = \frac{\lambda_e}{\lambda_e + \lambda_c}. \$\$

Given \\\lambda_e\\ and a target probability \\p\\, the function solves
for the censoring rate \\\lambda_c\\:

\$\$ \lambda_c = \frac{\lambda_e}{p} - \lambda_e. \$\$

In this mode, the returned value is the required independent censoring
rate.

## Mode = `"admin"`

Assume a single event time \\T \sim \mathrm{Exp}(\lambda_e)\\ with fixed
administrative censoring at time \\A \> 0\\.

The probability of observing the event by time \\A\\ is

\$\$ P(T \le A) = 1 - e^{-\lambda_e A}. \$\$

Given \\A\\ and a target probability \\p\\, the function solves for the
event rate \\\lambda_e\\:

\$\$ \lambda_e = -\frac{\log(1 - p)}{A}. \$\$

In this mode, the returned value is the required event rate.

## Mode = `"semi-competing"`

This mode uses a low-correlation or independent-clock approximation for
a secondary non-fatal event in a semi-competing risks setting.

Let:

- \\\lambda\_{e1}\\ be the fatal event rate,

- \\\lambda\_{c1}\\ be the censoring rate for the fatal endpoint,

- \\\lambda\_{e2}\\ be the non-fatal event rate, and

- \\\lambda\_{c2}\\ be the censoring rate for the non-fatal endpoint.

Under the approximation used in the package vignette, the probability
that the secondary non-fatal event is observed first is

\$\$ p \approx \frac{\lambda\_{e2}} {\lambda\_{e1} + \lambda\_{c1} +
\lambda\_{e2} + \lambda\_{c2}}. \$\$

Given \\p\\, \\\lambda\_{e1}\\, \\\lambda\_{c1}\\, and
\\\lambda\_{e2}\\, the function solves for \\\lambda\_{c2}\\:

\$\$ \lambda\_{c2} = \frac{\lambda\_{e2}}{p} - \left(\lambda\_{e1} +
\lambda\_{c1} + \lambda\_{e2}\right). \$\$

In this mode, the returned value is the approximate non-fatal censoring
rate.

Because this is an approximation, the realized observed non-fatal event
probability in
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
may differ once dependence, censoring, and fatal-event logic are fully
applied.

## See also

[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
for simulation of trial datasets using the resulting rates.

## Examples

``` r
# ------------------------------------------------------------
# 1) Independent exponential censoring only
# Solve for censoring rate given event rate and target event probability
# ------------------------------------------------------------
rate_from_prob(
  target_prob = 0.90,
  mode = "simple",
  event_rate = 1 / 24
)
#> [1] 0.00462963
# approximately 1/216

# ------------------------------------------------------------
# 2) Administrative censoring only
# Solve for event rate given follow-up time and target event probability
# ------------------------------------------------------------
rate_from_prob(
  target_prob = 0.20,
  mode = "admin",
  admin_time = 4
)
#> [1] 0.05578589
# approximately 0.0558

# ------------------------------------------------------------
# 3) Semi-competing risks approximation
# Solve for the non-fatal censoring rate
# ------------------------------------------------------------
rate_from_prob(
  target_prob = 0.20,
  mode = "semi-competing",
  fatal_event_rate = 1 / 50,
  fatal_censor_rate = 1 / 16.667,
  nonfatal_event_rate = 1 / 35
)
#> [1] 0.03428691
```
