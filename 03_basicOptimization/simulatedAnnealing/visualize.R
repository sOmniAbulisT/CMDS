# ==============================================================================
# Project: Simulated Annealing Optimization Engine
# File: visualize.R
# Description: Reads the saved SA results and plots the convergence trajectory
#              using ggplot2.
# ==============================================================================

library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Load the Simulated Annealing Results
# ------------------------------------------------------------------------------
resultFile <- "sa_result.rds"

if (!file.exists(resultFile)) {
  stop("Error: 'sa_result.rds' not found. Please run 'main.R' first to generate results.")
}

saResult <- readRDS(resultFile)

# ------------------------------------------------------------------------------
# 2. Prepare Data for Plotting
# Extract the history array. Remember C++ objPeaks uses a negative sign for 
# minimization, so we revert it back to the true maximum value.
# ------------------------------------------------------------------------------
historyMaxValues <- -saResult$history
iterations       <- 0:(length(historyMaxValues) - 1)

plotData <- data.frame(
  Iteration = iterations,
  MaxFitness = historyMaxValues
)

# ------------------------------------------------------------------------------
# 3. Create the Convergence Plot using ggplot2
# ------------------------------------------------------------------------------
cat("Generating convergence plot...\n")

convergencePlot <- ggplot(plotData, aes(x = Iteration, y = MaxFitness)) +
  geom_line(color = "#2c3e50", size = 1) +
  geom_point(data = tail(plotData, 1), aes(x = Iteration, y = MaxFitness), 
             color = "#e74c3c", size = 3) +
  labs(
    title = "Simulated Annealing Convergence Trajectory",
    subtitle = paste("Global Maximum Found at z =", round(tail(historyMaxValues, 1), 4)),
    x = "Temperature Step (Outer Loop Iteration)",
    y = "Objective Value (z)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(color = "#7f8c8d", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# 4. Save and Display the Plot
# ------------------------------------------------------------------------------
outputImage <- "sa_convergence.png"
ggsave(outputImage, plot = convergencePlot, width = 8, height = 5, dpi = 300)

cat(sprintf("Success! Convergence plot successfully saved as '%s'\n", outputImage))