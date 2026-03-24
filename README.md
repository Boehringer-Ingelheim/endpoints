Overview of `endpoints`
================

`endpoints` is an R package for flexibly simulating clinical trial
datasets with correlated endpoints of different types through a Gaussian
copula framework. Centered on a single core function, `makeData()`, it
provides a unified interface for generating continuous, binary, count,
and time-to-event outcomes, while also supporting multiple treatment
arms, censoring, and stochastic enrollment. The package is designed to
produce realistic patient-level datasets with configurable treatment
effects and dependence structures for trial planning and simulation
studies.

## Features

`endpoints` supports

- Multi-endpoint simulation of correlated clinical trial data within a
  unified framework
- Mixed endpoint types, including continuous, binary, count, and
  time-to-event outcomes
- Multiple treatment arms with arm-specific treatment effects
- Correlation control through a Gaussian copula, with optional automatic
  calibration to target observed Pearson correlations
- Time-to-event functionality including independent censoring,
  administrative censoring, and fatal/non-fatal endpoint logic for
  semi-competing risks settings
- Stochastic enrollment with uniform, exponential, and piecewise
  exponential accrual mechanisms
- Longitudinal data creation, either by treating repeated measurements
  as correlated endpoints or by generating correlated latent quantities
  for downstream longitudinal data-generating processes
- Comprehensive S3 methods for simulation objects, including print(),
  summary(), and plot()

## Installation

Install from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("Boehringer-Ingelheim/endpoints")
```

## Documentation

- [A User Guide for the `endpoints`
  Package](https://boehringer-ingelheim.github.io/endpoints/articles/user_guide.html)
- [Simulating Longitudinal Data with
  `endpoints`](https://boehringer-ingelheim.github.io/endpoints/articles/longitudinal_data.html)

## Quick Start

Generate a data set with three correlated endpoints, one time-to-event,
one continuous and one binary:

``` r
library(endpoints)

# correlation structure across the 3 endpoints
corMat_3 <- corr_make(
  num_endpoints = 3,
  values = rbind(
    c(1, 2, 0.20),
    c(1, 3, 0.10),
    c(2, 3, 0.15)
  )
)

# endpoint specifications
cont_ep <- list(
  endpoint_type = "continuous",
  baseline_mean = 10,
  sd            = c(3, 2),
  trt_effect    = -2
)

bin_ep <- list(
  endpoint_type = "binary",
  baseline_prob = 0.30,
  trt_prob      = 0.45
)

tte_ep <- list(
  endpoint_type  = "tte",
  baseline_rate  = 1 / 24,
  trt_effect     = log(0.80),
  censoring_rate = 1 / 216,
  fatal_event    = FALSE
)

# simulate data
sim_dat <- makeData(
  correlation_matrix    = corMat_3,
  sample_size_per_group = 1000,
  SEED                  = 123,
  endpoint_details      = list(cont_ep, bin_ep, tte_ep)
)

# inspect the simulation
sim_dat
summary(sim_dat)
head(sim_dat$data)

# quick pairwise plot for arm 0
plot(sim_dat, arm = 0, names = c("Biomarker", "Response", "Time to event"))
```
