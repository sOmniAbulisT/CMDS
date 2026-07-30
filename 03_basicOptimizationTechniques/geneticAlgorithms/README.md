# Genetic Algorithm Optimization

This directory implements a generalized, high-performance **Genetic Algorithm** engine built on a modular **R / C++ (Rcpp)** 
architecture. The solver leverages principles of natural selection and molecular genetics to efficiently navigate 
non-convex and discontinuous discrete landscapes.

---

## Algorithmic Architecture

Unlike local trajectory search methods (e.g., Simulated Annealing or Tabu Search) that track a single coordinate 
vector, this engine maintains and evolves an explicit collective **Population of Designs** ($m$) simultaneously. 

### Core Optimization Pillars:

1. **Modular Operator Design (`inline` Optimized):**
   The critical evolutionary steps—Selection, Crossover, and Mutation-are decoupled into independent, semantic inline functions in C++. This guarantees architectural clean code and maintainability while retaining structural $-O3$ register loop expansion speeds.
2. **Scale-Invariant Tournament Selection ($k=3$):**
   To bypass the classical absolute fitness translation/offset vulnerability common in Roulette-Wheel schemes, the mating pool is formed via Rank-based Tournament Selection. This keeps selection pressure robust and stable under arbitrary linear data offsets.
3. **Stochastic Genetic Recombination & Mutation:**
   * **Crossover:** Implements a single-point crossover mask activated by a baseline threshold probability ($p_c$).
   * **Mutation:** Operates via uniform bit-wise mutation driven by individual element-level probabilities ($p_m$), leveraging fast bitwise XOR primitives (`^= 1`) for structural bit-flipping.
4. **Dual Convergence Criteria (Statistical Guard):**
   The engine terminates gracefully when encountering either boundary condition:
   * **Iteration Limit:** The active generation marker reaches the absolute cap threshold ($i \ge i_{\text{max}}$).
   * **Population Trait Homogeneity:** The variance/standard deviation of fitness distributions across all active chromosomes drops below a strict stability margin ($s_f \le s_{f, \text{max}}$), indicating stable mathematical convergence.

---

## Directory Layout

```text
geneticAlgorithms/
├── README.md             # Technical manual and optimization engine specifications
├── main.R                # Master script: defines benchmark objectives, drives Rcpp, and saves RDS
├── visualize.R           # Post-processor: renders multi-generational fitness convergence profiles
└── src/
    ├── Makevars          # Compilation parameters (-O3 loop vectorization and C++17 flags)
    └── gaEngine.cpp      # Inline-optimized core evolutionary computing kernel
```

---

## Execution Routine

### Prerequisites
A operational modern C++ compiler toolchain (e.g., Rtools or GCC via WSL) must be mapped to your R execution 
environment, alongside the `Rcpp` and `ggplot2` visualization suites.

### Step 1: Initiate Evolutionary Synthesis
Source the primary routing console script to compile the C++ source files and execute the multi-generational search loop:

```R
source("03_basicOptimizationTechniques/geneticAlgorithms/main.R")
```
This generates the randomized initial population, triggers the C++ modular matrix mechanics, logs statistics to stdout, and serializes analytics 
into `ga_result.rds`.

### Step 2: Extract Post-Run Convergence Visualization
Trigger the graphics module to review generational performance and variance trajectories:

```R
source("03_basicOptimizationTechniques/geneticAlgorithms/visualize.R")
```
This processes the serialized array matrices and exports a high-resolution convergence path plot as 
`ga_convergence.png`.


