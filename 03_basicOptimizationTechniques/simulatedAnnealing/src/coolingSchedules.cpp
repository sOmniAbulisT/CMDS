#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;

//' Geometric (Exponential) cooling
//' 
// [[Rcpp::export]]
double tempGeometric (double initialTemp, double alpha, int iteration)
{
  return initialTemp * std::pow(alpha, iteration); 
}

//' Linear cooling
//' 
// [[Rcpp::export]]
double tempLinear (double initialTemp, double alpha, int iteration)
{
  double currentTemp = initialTemp - alpha * iteration;
  return (currentTemp > 0.0) ? currentTemp : 1e-7; 
}

//' Logarithmic cooling
//' 
// [[Rcpp::export]]
double tempLogarithmic (double initialTemp, double alpha, int iteration)
{
  return initialTemp / std::log(1.0 + iteration); 
}