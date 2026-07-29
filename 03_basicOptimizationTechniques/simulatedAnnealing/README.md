# Simulated Annealing Optimization

This repository contains a high-performance implementation of the **Simulated Annealing** algorithm. Designed with a modular architecture, the core computational engine is built in **C++ (via RcppArmadillo)** for maximum speed, while data processing, execution, and visualization are handled in **R**.

## Directory Overview
Simulated Annealing is a probabilistic technique for approximating the global optimum of a given function, heavily inspired by the physical process of metallurgy. This module is capable of escaping local optima to solve complex, highly non-convex, and multi-modal optimization problems.

**Key Features:**
* **C++ Core Engine:** Employs a double-loop architecture (Markov Chain inner loop + Cooling Schedule outer loop) for optimal performance.
* **Modular Design:** Objective functions, perturbation mechanisms, and temperature schedules are decoupled for easy extension to industrial scheduling or experimental design (DOE) problems.
* **Robustness:** Independent of initial guesses and functional continuity/convexity.

## Directory Structure
* `notes/`: Quarto (`.qmd`) documents containing mathematical formulations (Metropolis criterion, Boltzmann distribution) and algorithm theoretical background.
* `src/`: C++ source files containing the SA engine, cooling schedules, and target objective functions.
* `main.R`: Main entry point to configure parameters, compile the C++ core, and execute the solver.
* `visualize.R`: `ggplot2` scripts for plotting cooling curves and convergence trajectories.

## Quick Start
1. Ensure you have R and `RcppArmadillo` installed.
2. Clone this repository.
3. Source the main execution script in your R console:

```R
source("03_basicOptimizationTechniques/simulatedAnnealing/main.R")
```