#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//[[Rcpp::export]]
List PowerIteration (NumericMatrix A, int maxIter = 1000, double tol = 1e-12)
{
  int n = A.nrow();
  if (n != A.ncol()) stop("Matrix must be square");
  
  NumericVector v(n); 
  
  // Initialization
  for (int i = 0; i < n; i++)
  {
    v[i] = R::rnorm(0.0, 1.0); 
  }
  
  // L2 normalization
  double norm = 0.0; 
  for(int i = 0; i < n; i++) norm += v[i] * v[i];
  norm = std::sqrt(norm);
  for(int i = 0; i < n; i++) v[i] /= norm;
  
  double lambdaOld = 0.0; 
  double lambdaNew = 0.0; 
  
  // Main Loop
  for (int k = 0; k < maxIter; k++)
  {
    NumericVector w(n); 
    
    // matrix multiplier
    for (int i = 0; i < n; i++)
    {
      double sum = 0.0; 
      for (int j = 0; j < n; j++)
      {
        sum += A(i, j) * v[j]; 
      }
      w[i] = sum; 
    }
    
    // Rayleigh Quotient
    lambdaNew = 0.0; 
    for (int i = 0; i < n; i++)
    {
      lambdaNew += v[i] * w[i]; 
    }
    
    if (std::abs(lambdaNew - lambdaOld) < tol) 
    {
      break;
    }
    lambdaOld = lambdaNew;
    
    norm = 0.0;
    for(int i = 0; i < n; i++) norm += w[i] * w[i];
    norm = std::sqrt(norm);
    for(int i = 0; i < n; i++) v[i] = w[i] / norm;
  }
  
  return List::create(
    Named("eigenvalue") = lambdaNew,
    Named("eigenvector") = v
  );
}