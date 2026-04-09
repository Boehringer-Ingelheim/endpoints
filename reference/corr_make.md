# Construct a correlation matrix from endpoint-index triplets

`corr_make()` constructs a correlation matrix from a set of \\(i, j,
\rho)\\ triplets, where \\i\\ and \\j\\ are endpoint indices and
\\\rho\\ is their desired pairwise correlation.

## Usage

``` r
corr_make(num_endpoints, values = NULL)
```

## Arguments

- num_endpoints:

  A single positive integer giving the number of endpoints, that is, the
  dimension of the correlation matrix to construct.

- values:

  Optional object coercible to a numeric matrix with 3 columns. Each row
  must have the form `c(i, j, rho)`, where:

  i

  :   First endpoint index.

  j

  :   Second endpoint index.

  rho

  :   Desired correlation between endpoints `i` and `j`.

  `values` may be a matrix, data frame, or vector coercible to 3
  columns.

  The endpoint indices in `values` correspond to the order in which
  endpoints are supplied in `endpoint_details`; for example, index 1
  refers to the first endpoint in the list, index 2 to the second, and
  so on.

  If `NULL`, the identity matrix of size `num_endpoints` is returned.
  Unspecified off-diagonal entries default to 0.

## Value

A symmetric numeric correlation matrix of dimension
`num_endpoints x num_endpoints`, with ones on the diagonal.

## Details

This is a convenience function for specifying sparse correlation
structures when using
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md).
Any unspecified off-diagonal entries are left at 0, and diagonal entries
are set to 1.

## How entries are interpreted

For each row in `values`, the indices `i` and `j` are mapped to the
upper triangle of the matrix using: \$\$ i^\star = \min(i, j), \qquad
j^\star = \max(i, j), \$\$ so the order of the indices does not matter.
The specified correlation is then assigned symmetrically: \$\$
R\[i^\star, j^\star\] = R\[j^\star, i^\star\] = \rho. \$\$ Finally, the
diagonal is reset to exactly 1, so any user-supplied diagonal entries
are overwritten.

## Use with [`makeData()`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)

The resulting matrix can be passed directly to
[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
as `correlation_matrix`.

When `target_correlation = FALSE`, the resulting matrix is interpreted
as the latent Gaussian correlation matrix.

When `target_correlation = TRUE`, the resulting matrix is treated as the
desired observed Pearson correlation matrix, and the simulation routine
attempts to calibrate a latent Gaussian matrix to approximately match it
after marginal transformation.

## See also

[`makeData`](https://boehringer-ingelheim.github.io/endpoints/reference/makeData.md)
for the main simulation function.

[`endpoint_details`](https://boehringer-ingelheim.github.io/endpoints/reference/endpoint_details.md)
for endpoint specification details.

[`calibration_control`](https://boehringer-ingelheim.github.io/endpoints/reference/calibration_control.md)
for correlation calibration settings.

## Examples

``` r
## Identity correlation matrix
corr_make(num_endpoints = 3)
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1

## Three endpoints with pairwise specification
R3 <- corr_make(
  num_endpoints = 3,
  values = rbind(
    c(1, 2, 0.20),
    c(1, 3, 0.10),
    c(2, 3, 0.15)
  )
)
R3
#>      [,1] [,2] [,3]
#> [1,]  1.0 0.20 0.10
#> [2,]  0.2 1.00 0.15
#> [3,]  0.1 0.15 1.00

## AR(1)-style structure across five repeated measures
rho <- 0.5
R5 <- corr_make(
  num_endpoints = 5,
  values = rbind(
    c(1, 2, rho),   c(1, 3, rho^2), c(1, 4, rho^3), c(1, 5, rho^4),
    c(2, 3, rho),   c(2, 4, rho^2), c(2, 5, rho^3),
    c(3, 4, rho),   c(3, 5, rho^2),
    c(4, 5, rho)
  )
)
R5
#>        [,1]  [,2] [,3]  [,4]   [,5]
#> [1,] 1.0000 0.500 0.25 0.125 0.0625
#> [2,] 0.5000 1.000 0.50 0.250 0.1250
#> [3,] 0.2500 0.500 1.00 0.500 0.2500
#> [4,] 0.1250 0.250 0.50 1.000 0.5000
#> [5,] 0.0625 0.125 0.25 0.500 1.0000

## Order of indices does not matter
corr_make(
  num_endpoints = 3,
  values = rbind(
    c(2, 1, 0.30),
    c(3, 1, 0.10)
  )
)
#>      [,1] [,2] [,3]
#> [1,]  1.0  0.3  0.1
#> [2,]  0.3  1.0  0.0
#> [3,]  0.1  0.0  1.0
```
