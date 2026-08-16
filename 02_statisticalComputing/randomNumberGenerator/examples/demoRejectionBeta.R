## loading functions
source("02_statisticalComputing/randomNumberGenerator/R/generateRV.R")

# target pdf: beta(2, 2)
# f(x) = (gamma(4) / (gamma(2) * gamma(2))) * x^(2-1) * (1 - x)^(2-1) = 6*x*(1-x)
f <- function(x){
  return(6 * x * (1 - x))
}

# proposal pdf g(x): uniform(0, 1)
# Since the value range of Beta is between 0 and 1, uniform distribution is the most suitable base.
g <- function(x){
  return(dunif(x, min = 0, max = 1))
}

# proposal generator
rg <- function(n){
  return(runif(n, min = 0, max = 1))
}

# bounding constant
# c * g(x) >= f(x)
c <- 1.5

betaSample <- GenerateRV(n = 10000, method = "rejection", f = f, g = g, rg = rg, c = c)

## ploting
png("02_statisticalComputing/randomNumberGenerator/examples/plots/rejectionBeta.png", width = 800, height = 600, res = 120)

hist(betaSample, breaks = 50, probability = TRUE, 
     main = "Rejection Sampling: Beta(2, 2)", 
     xlab = "x", col = "lightblue", border = "white")

curve(dbeta(x, 2, 2), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = c("Empirical Density", "Theoretical PDF"), 
       fill = c("lightblue", NA), border = c("black", NA), 
       col = c(NA, "red"), lwd = c(NA, 2))

dev.off()
