## loading functions
source("02_statisticalComputing/randomNumberGenerator/R/generateRV.R")

# target pdf: Exponential(rate = 2)
# inverse cdf: F^{-1}(u) = -log(1 - u) / 2
invCDFexp <- function(u){
  return(-log(1 - u) / 2)
}

##
expSample <- GenerateRV(n = 10000, method = "inverse", InvCDF = invCDFexp)

## ploting
png("02_statisticalComputing/randomNumberGenerator/examples/plots/inverseExp.png", width = 800, height = 600, res = 120)

hist(expSample, breaks = 50, probability = TRUE, 
     main = "Inverse Transform: Exponential(rate = 2)", 
     xlab = "x", col = "lightblue", border = "white")

curve(dexp(x, rate = 2), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = c("Emperical Density", "Theoratical PDF"), 
       fill = c("lightblue", NA), border = c("black", NA), 
       col = c(NA, "red"), lwd = c(NA, 2))

dev.off()
