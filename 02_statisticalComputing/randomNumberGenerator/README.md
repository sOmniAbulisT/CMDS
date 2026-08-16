# Universal Random Variable Generator

## Overview
This repository contains a foundational statistical computing project focused on **Stochastic Simulation and Sampling Methods**. It features a custom-built, universal R function (`GenerateRV`) that generates random variables from scratch without relying on high-level, built-in distribution functions. 

By manually implementing core algorithmic logic, this project demonstrates a deep understanding of probability distributions, algorithmic optimization, and R programming vectorization—essential skills for advanced statistical modeling and computational research.

## Core Function: `GenerateRV`
The heart of this project is the `GenerateRV()` wrapper function. It acts as a unified interface to draw random samples using five classical statistical algorithms.

### Supported Algorithms
*   **Inverse Transform Method**: Generates samples by computing the inverse of the target Cumulative Distribution Function (CDF).
*   **Acceptance-Rejection Method**: Utilizes a proposal distribution and a bounding constant to sample from complex Probability Density Functions (PDFs).
*   **Transformation Method**: Applies mathematical transformations to base random variables.
*   **Sums Method**: Generates variables based on the sum of underlying distributions (e.g., Central Limit Theorem applications).
*   **Mixture Method**: Employs a two-stage sampling process to generate data from finite mixture models.

## File Structure

```
randomNumberGenerator/
├── R/
│   └── GenerateRV.R              # Core wrapper function and algorithm logic
├── examples/
│   ├── demo_inverse_exp.R        # Demo: Exponential distribution via Inverse Transform
│   ├── demo_rejection_beta.R     # Demo: Beta distribution via Acceptance-Rejection
│   ├── demo_mixture_dist.R       # Demo: Gaussian Mixture distribution
│   └── plots/                    # Visualizations of generated distributions
└── README.md                     
```

## Usage and Examples
Here is a demonstration of using `GenerateRV` to sample from a Gaussian Mixture Distribution 
(30% $N(-3, 1)$ and 70% $N(3, 1)$).

```{R}
# Source the core function
source("R/GenerateRV.R")

# Define Stage-1: Population parameter generator
stage1_rg <- function(n) {
  sample(c(-3, 3), size = n, replace = TRUE, prob = c(0.3, 0.7))
}

# Define Stage-2: Conditional normal generator
stage2_rg <- function(n, params) {
  rnorm(n, mean = params, sd = 1)
}

# Generate 10,000 samples
mix_samples <- GenerateRV(n = 10000, method = "mixture", rg1 = stage1_rg, rg2 = stage2_rg)
```

## Visualization
The empirical density of the generated samples perfectly matches the theoretical probability density function:

(Note: Ensure the plot image is located at `examples/plots/mixture_bimodal.png`)