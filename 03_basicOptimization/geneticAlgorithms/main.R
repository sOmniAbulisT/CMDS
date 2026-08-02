library(Rcpp)

cat(" Compiling C++ GA Engine...\n")
Rcpp::sourceCpp("03_basicOptimizationTechniques/geneticAlgorithms/src/gaEngine.cpp")

decode_binary <- function(genes, bits_per_var, min_val, max_val) {
  n_vars <- length(genes) / bits_per_var
  decoded_vals <- numeric(n_vars)
  
  for (i in 1:n_vars) {
    start_idx <- (i - 1) * bits_per_var + 1
    end_idx <- i * bits_per_var
    var_bits <- genes[start_idx:end_idx]
    
    decimal_val <- sum(var_bits * (2^(seq_along(var_bits) - 1)))
    max_decimal <- (2^bits_per_var) - 1
    decoded_vals[i] <- min_val + (decimal_val / max_decimal) * (max_val - min_val)
  }
  return(decoded_vals)
}

rastrigin_objective <- function(genes) {
  x <- decode_binary(genes, bits_per_var = 10, min_val = -5.12, max_val = 5.12)
  A <- 10
  fitness <- A * length(x) + sum(x^2 - A * cos(2 * pi * x))
  return(fitness)
}

POP_SIZE        <- 100
BITS_PER_VAR    <- 10
N_VARS          <- 2
CHROM_LENGTH    <- N_VARS * BITS_PER_VAR
P_CROSSOVER     <- 0.8
P_MUTATION      <- 0.05
MAX_GENERATIONS <- 200
MIN_STD_DEV     <- 1e-4

cat("Running GA Solver...\n")
set.seed(42)

start_time <- Sys.time()
ga_result <- gaEngine(
  objFunc        = rastrigin_objective,
  popSize        = POP_SIZE,
  chromLength    = CHROM_LENGTH,
  pCrossover     = P_CROSSOVER,
  pMutation      = P_MUTATION,
  maxGenerations = MAX_GENERATIONS,
  minStdDev      = MIN_STD_DEV
)
end_time <- Sys.time()

final_x <- decode_binary(ga_result$bestGenes, bits_per_var = BITS_PER_VAR, min_val = -5.12, max_val = 5.12)

cat(sprintf("\n Done! Runtime: %.4fs | Gens: %d\n", as.numeric(end_time - start_time, units="secs"), ga_result$generationsRun))
cat(sprintf("   Best Fitness: %.6f\n", ga_result$bestFitness))
cat(sprintf("   Coordinates : [%.4f, %.4f]\n\n", final_x[1], final_x[2]))

saveRDS(ga_result, file = "ga_result.rds")