# Specification format for `calibration_control` used by `makeData()`

`calibration_control` is a named list of tuning parameters used by
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
when `target_correlation = TRUE`.

## Details

These settings control the numerical search used to calibrate a latent
Gaussian correlation matrix so that, after transformation through the
endpoint-specific marginals, the observed Pearson correlation matrix is
approximately equal to the user-supplied `correlation_matrix`.

## Default values

The default list is:

     list( n_mc =
    10000, tol = 0.001, maxit = 100, rho_cap = 0.999, ensure_pd = TRUE,
    conv_norm_type = "F" ) 

## Arguments

- list("n_mc"):

  Monte Carlo sample size used internally when approximating the
  transformed correlation induced by a candidate latent Gaussian
  correlation. Larger values generally improve numerical stability but
  increase computation time.

- list("tol"):

  Positive numeric tolerance passed to the one-dimensional root-finding
  routine used in pairwise calibration. Smaller values can improve
  accuracy but may require more iterations and computation.

- list("maxit"):

  Positive integer giving the maximum number of root-finding iterations.

- list("rho_cap"):

  Numeric scalar in (0,1) giving the maximum absolute latent Gaussian
  correlation considered during calibration. This avoids numerical
  instability at the exact boundaries \\(\pm 1\\).

- list("ensure_pd"):

  Logical scalar. If `TRUE`, the calibrated latent correlation matrix is
  repaired (if needed) to ensure it is a valid positive-definite
  correlation matrix before simulation proceeds.

- list("conv_norm_type"):

  Character scalar passed through to positive-definiteness repair
  routines (for example, `Matrix::nearPD(..., conv.norm.type = ...)`
  when available). The default `"F"` (Frobenius)).

## When calibration is used

These settings are used only when:

- `correlation_matrix` is not `NULL`, and

- `target_correlation = TRUE`.

If `target_correlation = FALSE`, the supplied `correlation_matrix` is
used directly as the latent Gaussian correlation matrix, and
`calibration_control` is ignored. In general, it is recommended that
`target_correlation = TRUE`.

## Practical guidance

- In many applications, the defaults have been observed to be
  sufficient.

- If the estimated correlation matrices differ substantially from the
  target, the increasing `n_mc` will usually have the most impact. Note
  that increasing this value will increase computation time.

- When all simulated endpoints are continuous Gaussian,
  target_correlation = FALSE may be used as Pearson correlation is
  preserved under the normal marginal transformation and no calibration
  is needed. This may be especially important in simulation studies as
  this will save computation time. %

- For publication-quality simulation studies where observed Pearson
  correlation matching is important, % increasing `n_mc` may improve
  stability. %

- If calibration becomes slow, consider reducing `n_mc`, loosening
  `tol`, or using % `target_correlation = FALSE` when the latent
  Gaussian correlation is acceptable.

- Even with calibration, finite-sample summaries from the simulated data
  will not match the target matrix exactly.

- Please note that not all combinations of endpoint marginals and target
  Pearson correlations are feasible, so some requested correlation
  structures may be unattainable under the specified data-generating
  distributions.

## See also

[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
for the main simulation routine.

[`endpoint_details`](https://boehringer-ingelheim.github.io/endpoints/reference/endpoint_details.md)
for endpoint definitions.

[`enrollment_details`](https://boehringer-ingelheim.github.io/endpoints/reference/enrollment_details.md)
for follow-up and enrollment settings.

## Examples

``` r
## Default settings
calibration_control <- list(
  n_mc = 10000,
  tol = 0.001,
  maxit = 100,
  rho_cap = 0.999,
  ensure_pd = TRUE,
  conv_norm_type = "F"
)

## Higher-precision example
calibration_control <- list(
  n_mc = 50000,
  tol = 1e-4,
  maxit = 200,
  rho_cap = 0.999,
  ensure_pd = TRUE,
  conv_norm_type = "F"
)
```
