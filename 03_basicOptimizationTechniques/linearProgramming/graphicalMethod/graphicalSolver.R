# ==============================================================================
# Project: Linear Programming Solver (Graphical Method)
# File: graphicalSolver.R
# Author: Kun-Hong Liao
# Description: 2D visualization of the feasible region and optimal solution
#              for the Resource Allocation Problem using ggplot2.
# ==============================================================================

# Load required library (install if necessary: install.packages("ggplot2"))
library(ggplot2)

cat("Initializing graphical solver...\n")

# 1. Define the vertices of the feasible region
# Mathematical calculation of intersection points for the binding constraints:
# P1: Origin (0, 0)
# P2: Intersection of C3 (Demand) and x-axis -> (15, 0)
# P3: Intersection of C1 (Raw Material) and C3 -> x1 = 15, x2 = (21-15)/4 = 1.5
# P4: Intersection of C1 and C2 (Capacity) -> x1 = 63/11 (~5.73), x2 = 42/11 (~3.82)

vertices <- data.frame(
  x1 = c(0, 15, 15, 63/11),
  x2 = c(0, 0, 1.5, 42/11),
  label = c("P1 (0, 0)", "P2 (15, 0)", "P3 (15, 1.5)", "P4 (5.73, 3.82)")
)

# Calculate Objective Function Z = 50*x1 + 200*x2 at each vertex
vertices$Z <- 50 * vertices$x1 + 200 * vertices$x2

# Identify the optimal vertex (Highest Z)
optimalPoint <- vertices[which.max(vertices$Z), ]

# 2. Prepare data for plotting the constraint lines
xVal <- seq(0, 20, length.out = 100)

# C1: x1 + 4x2 = 21  => x2 = (21 - x1) / 4
# C2: -4x1 + 6x2 = 0 => x2 = (4/6) * x1
linesData <- data.frame(
  x1 = rep(xVal, 2),
  x2 = c((21 - xVal) / 4, (4/6) * xVal),
  Constraint = rep(c("C1: x1 + 4x2 = 21 (Raw Material)", 
                     "C2: -4x1 + 6x2 = 0 (Capacity)"), each = length(xVal))
)

# Filter out negative y-values for a cleaner plot
linesData <- subset(linesData, x2 >= 0)

# 3. Generate the ggplot2 visualization
cat("Rendering the feasible region plot...\n")

p <- ggplot() +
  # Draw the Feasible Region polygon
  geom_polygon(data = vertices, aes(x = x1, y = x2), fill = "steelblue", alpha = 0.3) +
  
  # Draw Constraint Lines
  geom_line(data = linesData, aes(x = x1, y = x2, color = Constraint), linewidth = 1) +
  geom_vline(aes(xintercept = 15, color = "C3: x1 = 15 (Demand)"), linewidth = 1, linetype = "dashed") +
  
  # Highlight the vertices
  geom_point(data = vertices, aes(x = x1, y = x2), size = 3, color = "black") +
  geom_text(data = vertices, aes(x = x1, y = x2, label = label), 
            hjust = -0.15, vjust = -0.5, size = 4) +
  
  # Highlight the Optimal Solution
  geom_point(data = optimalPoint, aes(x = x1, y = x2), size = 5, color = "firebrick") +
  geom_text(data = optimalPoint, aes(x = x1, y = x2, label = paste0("Optimal Z = ", round(Z, 2))),
            hjust = 1.1, vjust = -1, color = "firebrick", fontface = "bold", size = 5) +
  
  # Axes and Labels customization
  scale_x_continuous(limits = c(0, 20), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 8), expand = c(0, 0)) +
  labs(
    title = "Linear Programming: Graphical Method",
    subtitle = paste0("Optimal Solution found at x1 = ", round(optimalPoint$x1, 2), 
                      ", x2 = ", round(optimalPoint$x2, 2), " | Max Revenue Z = ", round(optimalPoint$Z, 2)),
    x = bquote(x[1] ~ " (5mm Cable in thousands of meters)"),
    y = bquote(x[2] ~ " (10mm Cable in thousands of meters)"),
    color = "Constraint Boundaries"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  )

# Display the plot
print(p)
cat("Done! Check your Plots pane.\n")
