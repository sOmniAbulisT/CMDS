# ==============================================================================
# Project: Simulated Annealing Optimization Engine
# File: main.R
# Description: Main entry point for compiling C++ core, setting hyperparameters,
#              and executing Simulated Annealing on the MIT Peaks Function.
# ==============================================================================

library(Rcpp)
library(RcppArmadillo)

# ------------------------------------------------------------------------------
# 1. Compile C++ Core Modules
# ------------------------------------------------------------------------------
cat("Compiling C++ core modules...\n")
sourceCpp("03_basicOptimizationTechniques/simulatedAnnealing/src/coolingSchedules.cpp")
sourceCpp("03_basicOptimizationTechniques/simulatedAnnealing/src/objectiveFuncs.cpp")
sourceCpp("03_basicOptimizationTechniques/simulatedAnnealing/src/saEngine.cpp")
cat("C++ modules successfully compiled!\n\n")

# ------------------------------------------------------------------------------
# 2. Hyper-parameter Configuration
# Based on the MIT lecture example: starting at (-2, -2) to find the global 
# maximum of the Peaks function.
# ------------------------------------------------------------------------------
initState   <- c(-2.0, -2.0)  # Initial search starting point
initialTemp <- 100.0          # Initial high temperature (enables early random walk)
alpha       <- 0.95           # Geometric cooling rate
minTemp     <- 1e-4           # Minimum temperature termination criterion
markovIter  <- 50             # Markov chain length (iterations to reach thermal equilibrium)

# ------------------------------------------------------------------------------
# 3. Execute Simulated Annealing Algorithm
# ------------------------------------------------------------------------------
cat("==================================================\n")
cat("Executing Simulated Annealing algorithm...\n")
cat("==================================================\n")

startTime <- Sys.time()

# Call the high-performance C++ engine
saResult <- SA(
  initState   = initState,
  objFunc     = objPeaks,
  perturbFunc = perturbNormal,
  initialTemp = initialTemp,
  alpha       = alpha,
  minTemp     = minTemp,
  markovIter  = markovIter
)

endTime <- Sys.time()
executionTime <- difftime(endTime, startTime, units = "secs")

# ------------------------------------------------------------------------------
# 4. Parse and Output Results
# Since C++ objPeaks converts the maximization problem to minimization 
# (via a negative sign), the true maximum = -bestEnergy
# ------------------------------------------------------------------------------
foundState  <- saResult$bestState
foundMaxVal <- -saResult$bestEnergy

cat(sprintf("Execution Time: %.4f seconds\n", executionTime))
cat(sprintf("Initial Starting Point: (%.3f, %.3f)\n", initState[1], initState[2]))
cat(sprintf("Global Maximum Located at (x, y): (%.4f, %.4f)\n", foundState[1], foundState[2]))
cat(sprintf("Global Maximum Value (z): %.4f\n", foundMaxVal))
cat("==================================================\n")

# ------------------------------------------------------------------------------
# 5. Save Results for Convergence Visualization
# ------------------------------------------------------------------------------
saveRDS(saResult, file = "sa_result.rds")
cat("Convergence results successfully saved to 'sa_result.rds'\n")