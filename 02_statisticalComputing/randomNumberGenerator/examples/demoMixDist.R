## loading functions
source("02_statisticalComputing/randomNumberGenerator/R/generateRV.R")

# target pdf: mixed normal distribution (30% N(3, 1), 70% N(3, 1))

# Stage 1: generating population parameter
stage1 <- function(n){
  mu <- sample(c(-3, 3), size = n, replace = TRUE, prob = c(0.3, 0.7))
  return(mu)
}

# Stage 2: Given the mean parameter, generate the corresponding normal random number.
stage2 <- function(n, params){
  return(rnorm(n, mean = params, sd = 1))
}

##
mixSample <- GenerateRV(n = 10000, method = "mixture", rg1 = stage1, rg2 = stage2)

## ploting
png("02_statisticalComputing/randomNumberGenerator/examples/plots/mixDist.png", width = 800, height = 600, res = 120)

hist(mixSample, breaks = 50, probability = TRUE, 
     main = "Mixture Distribution: Gaussian Mixture", 
     xlab = "x", col = "plum", border = "white")

theoreticalMix <- function(x) {
  0.3 * dnorm(x, mean = -3, sd = 1) + 0.7 * dnorm(x, mean = 3, sd = 1)
}

curve(theoreticalMix(x), add = TRUE, col = "red", lwd = 2)
legend("topleft", legend = c("Empirical Density", "Theoretical PDF"), 
       fill = c("plum", NA), border = c("black", NA), 
       col = c(NA, "red"), lwd = c(NA, 2))

dev.off()
