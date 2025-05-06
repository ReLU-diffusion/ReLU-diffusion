# ReLU Diffusion Model for "Go/No-Go" Decisions

This repository contains the R code for simulating and evaluating the ReLU diffusion model proposed in the NeurIPS 2025 submission titled *"A ReLU Diffusion Model for 'Go/No-Go' Decisions."*

## Files

- `simu_main.R`: Simulates repeated trials under the ReLU diffusion model, fits the model using method-of-moments (MOM), and reports the parameter estimation performance over multiple replications. It uses the `nleqslv` package for solving nonlinear equations.

- `Laplace_Transformation.R`: Computes the Laplace-transformed density and cumulative distribution function (CDF) of the model's stopping time using inverse Laplace techniques (via the `pracma` package) for various values of the drift parameter μ.

## Requirements

- R (version ≥ 4.0)
- R packages: `nleqslv`, `pracma`

To install required packages:

```r
install.packages("nleqslv")
install.packages("pracma")
```

## Running the Simulation

To reproduce the simulation results in the paper (Table 1), run:

```r
source("simu_main.R")
```

This will print out the mean, standard deviation, and RMSE for estimated parameters over 100 repetitions.

## Reproducing the Theoretical Distributions

To generate plots of the CDF or PDF under different drift values, run:

```r
source("Laplace_Transformation.R")
```

This script plots the CDF and PDFs under μ = 0, 0.5, and 1, replicating Figure 1 in the paper.

## License

MIT License
