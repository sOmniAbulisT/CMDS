---

editor_options: 
  markdown: 
    wrap: 72
---

# Linear Programming

Linear programming is a foundational mathematical method used to determine the best possible outcome (such as maximum profit or lowest cost) in a given mathematical model where the requirements are represented by linear relationships. It is widely used in resource allocation, operations research, and as a building block for more complex statistical and machine learning algorithms.

## Directory Structure

This module is structured to separate mathematical theory from high-performance algorithmic implementation:

``` text
linearProgramming/
├── README.md                     # Module overview (this file)
├── notes/                        
│   └── theoryFormulation.qmd     # Theory notes: formulation, standard form, and examples
├── graphicalMethod/
│   └── graphicalSolver.R         # R script for 2D visualization of feasible regions and vertices
├── src/                          
│   └── optimize.cpp              # Core C++ implementation using RcppArmadillo (Simplex Method)
└── main.R                        # Entry point: constructs the matrix and calls the C++ solver
```

## Mathematical Formulation (Standard Form)

To solve a linear programming problem computationally, real-world constraints must be translated into standard matrix notation.

The standard form of an linear programming problem is define as:

Objective function: $ \min (\text{or } \max) \quad Z = \mathbf{c}^T \mathbf{x} $

Subject to Constraints: $\mathbf{A}\mathbf{x} \le \mathbf{b}$ $\mathbf{x} \ge \mathbf{0}$

Where: \* $\mathbf{x}$ is the Design Vector (decision variables). \* $\mathbf{c}$ is the cost/weight vector. \* $\mathbf{A}$ is the constraint coefficient matrix. \* $\mathbf{b}$ is the resource bound vector.

By introducing Slack Variables, inequality constraints are converted into equations to form the initial Simplex Tableau, which is then solved iteratively using the Simplex Method.

## How to Run

### Prerequisites

Make sure you have R installed along with the following packages: \* `Rcpp` \* `RcppArmadillo`

### Execution

1.  Open the project in your preferred IDE (e.g., VS Code or RStudio).

2.  Run the `main.R` script.

3.  The script will automatically compile `src/optimize.cpp` and output the step-by-step matrix transformations and the final optimal solution.

This module is part of the graduate-level statistics and optimization coursework.
