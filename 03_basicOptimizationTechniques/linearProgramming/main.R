# ==============================================================================
# Project: Linear Programming Solver (Simplex Method)
# File: main.R
# Description: Entry point script to solve Linear Programming problems by
#              interfacing R with the C++ Simplex core solver via Rcpp.
# ==============================================================================

library(Rcpp)
sourceCpp("03_basicOptimizationTechniques/linearProgramming/src/optimize.cpp")

initialTableau <- matrix(c(
  1,   4,  1,  0,  0,  0,  21,   # Constraint 1: Raw material limit (x1 + 4x2 + s1 = 21)
  -4,   6,  0,  1,  0,  0,   0,   # Constraint 2: Capacity ratio (-4x1 + 6x2 + s2 = 0)
  1,   0,  0,  0,  1,  0,  15,   # Constraint 3: Market demand for 5mm (x1 + s3 = 15)
  -50, -200, 0,  0,  0,  1,   0    # Objective Row: Max Z = 50x1 + 200x2 -> -50x1 - 200x2 + Z = 0
), nrow = 4, byrow = TRUE)

colnames(initialTableau) <- c("x1", "x2", "s1", "s2", "s3", "Z", "RHS")
rownames(initialTableau) <- c("RawMaterial", "Capacity", "Demand", "Objective")

result <- runSimplexSolver(initialTableau)
