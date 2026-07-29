#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' 1. MIT Peaks Function
// [[Rcpp::export]]
double objPeaks(NumericVector currentState) 
{
   double xVal = currentState[0];
   double yVal = currentState[1];
   
   double termOne = 3.0 * std::pow(1.0 - xVal, 2) * std::exp(-std::pow(xVal, 2) - std::pow(yVal + 1.0, 2));
   double termTwo = -10.0 * (xVal / 5.0 - std::pow(xVal, 3) - std::pow(yVal, 5)) * std::exp(-std::pow(xVal, 2) - std::pow(yVal, 2));
   double termThree = -1.0 / 3.0 * std::exp(-std::pow(xVal + 1.0, 2) - std::pow(yVal, 2));
   
   double zVal = termOne + termTwo + termThree;
   
   return -zVal;
}
 
//' 2. Gaussian Perturbation
// [[Rcpp::export]]
NumericVector perturbNormal(NumericVector currentState, double stepSize = 0.5) 
{
   NumericVector newState = clone(currentState);
   
   for(int i = 0; i < newState.size(); ++i) 
   {
      newState[i] += R::rnorm(0.0, stepSize);
   }
   
   return newState;
}