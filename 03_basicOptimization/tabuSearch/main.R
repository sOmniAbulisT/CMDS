# ==============================================================================
# Project: Generic Tabu Search Optimization Engine
# File: main.R
# Description: Master script that compiles the C++ engine kernel, defines a 
#              benchmark objective function, runs the Tabu Search solver, 
#              and serializes the results into RDS format.
# ==============================================================================

library(Rcpp)

# ------------------------------------------------------------------------------
# 1. Compile and Source C++ Engine Kernel
# ------------------------------------------------------------------------------
sourceCpp("03_basicOptimizationTechniques/tabuSearch/src/tabuSearchEngin.cpp")

# ------------------------------------------------------------------------------
# 2. Define Benchmark Objective Function & Hyper-parameters
# ------------------------------------------------------------------------------

vectorLength <- 10
initState    <- sample(1:vectorLength)

benchmarkObjFunc <- function(vectorState) {
  weights <- 1:length(vectorState)
  targetPattern <- length(vectorState):1
  return(sum(weights * (vectorState - targetPattern)^2))
}

# hyper-parameter setting
maxIterations <- 100
tabuTenure    <- 7  

# ------------------------------------------------------------------------------
# 3. Execute High-Performance Tabu Search Engine
# ------------------------------------------------------------------------------
cat("Executing Tabu Search Optimization Loop...\n")

executionTime <- system.time({
  tsResult <- tsEngine(
    initState     = initState,
    objFunc       = benchmarkObjFunc,
    maxIterations = maxIterations,
    tabuTenure    = tabuTenure
  )
})

# ------------------------------------------------------------------------------
# 4. Display Execution Metrics
# ------------------------------------------------------------------------------
cat("\n==================================================\n")
cat("          Tabu Search Execution Summary           \n")
cat("==================================================\n")
cat(sprintf("Elapsed Time      : %.4f seconds\n", executionTime["elapsed"]))
cat(sprintf("Initial Energy    : %.4f\n", benchmarkObjFunc(initState)))
cat(sprintf("Best Energy Found : %.4f\n", tsResult$bestEnergy))
cat("Initial Vector    : ", paste(initState, collapse = ", "), "\n")
cat("Best Vector Found : ", paste(tsResult$bestPoint, collapse = ", "), "\n")
cat("Total History Steps: ", length(tsResult$history), "\n")
cat("==================================================\n\n")

# ------------------------------------------------------------------------------
# 5. Serialize Results to RDS
# ------------------------------------------------------------------------------
resultFile <- "ts_result.rds"
saveRDS(tsResult, file = resultFile)
cat(sprintf("Success! Results serialized and saved to '%s'\n", resultFile))