#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <numeric>   
#include <algorithm> 

using namespace Rcpp; 

// QR algorithm for eigenvalue solve and matrix multiply
#include "01_linearAlgebraicTechniques/eigenvalueProblem/qrEigenSolver.cpp"

// ==============================================================================
// Math Utilities
// ==============================================================================

// non-linear function
inline double gFunc (double u)
{
  return std::tanh(u); 
}

// derivative 
inline double gPrime (double u)
{
  double t = std::tanh(u); 
  return 1.0 - t * t; 
}

// normalization
void normalizeVector (NumericVector& v)
{
  double norm = 0.0; 
  for (int i = 0; i < v.size(); i++) norm += v[i] * v[i]; 
  norm = std::sqrt(norm); 
  for (int i = 0; i < v.size(); i++) v[i] /= (norm + 1e-12); 
}

// ==============================================================================
// Independent Component Analysis
// ==============================================================================
//[[Rcpp::export]]
List ICA (NumericMatrix X, int k, int maxIter = 200, double tol = 1e-5)
{
  int n = X.nrow();
  int p = X.ncol();
  
  if (k > p) stop("k must be less than p");
  
  // ---------------------------------------
  // Step 1: Centering
  // ---------------------------------------
  NumericMatrix Xc = clone(X); 
  for (int j = 0; j < p; j++)
  {
    double mean = 0.0; 
    for (int i = 0; i < n; i++) mean += Xc(i, j); 
    mean /= n;
    for (int i = 0; i < n; i++) Xc(i, j) -= mean; 
  }
  
  // ---------------------------------------
  // Step 2: Decorrelation
  // ---------------------------------------
  NumericMatrix Cov(p, p); 
  for (int i = 0; i < p; i++)
  {
    for (int j = i; j < p; j++)
    {
      double dot = 0.0; 
      for (int r = 0; r < n; r++) dot += Xc(r, i) * Xc(r, j); 
      Cov(i, j) = dot / (n - 1); 
      Cov(j, i) = Cov(i, j); 
    }
  }
  
  List eigenRes = QREigenSolver(Cov);
  NumericVector eigenvaluesRaw = eigenRes["eigenvalues"];
  NumericMatrix eigenvectorsRaw = eigenRes["eigenvectors"];
  
  std::vector<int> idx(p);
  std::iota(idx.begin(), idx.end(), 0);
  
  std::sort(idx.begin(), idx.end(), [&](int i, int j) 
  {
    return eigenvaluesRaw[i] > eigenvaluesRaw[j];
  });
  
  // -----------------------------------------
  // Step 3: Whitening 
  // -----------------------------------------
  NumericMatrix Vk(p, k); 
  NumerucVector dInvSqrt(k); 
  
  for (int j = 0; j < k; j++)
  {
    int sortedIdx = idx[j]; 
    
    dInvSqrt[j] = 1.0 / std::sqrt(eigenvaluesRaw[sortedIdx] + 1e-12); 
    
    for (int i = 0; i < p; i++)
    {
      Vk(i, k) = eigenvectorsRaw(i, sortedIdx) * dInvSqrt[j]; 
    }
  }
  
  NumericMatrix Z = matrixMultiply(Xc, Vk); 
  
  // ---------------------------------------
  // Step 4: FastICA iteration
  // ---------------------------------------
  NumericMatrix W(k, k); 
  
  for (int c = 0; c < k; c++)
  {
    NumericVector w(k); 
    for (int i = 0; i < k; i++) w[i] = R::rnorm(0.0, 1.0); 
    normalizeVector(w); 
    
    bool converged= false; 
    int iter = 0; 
    
    while (!converged && iter < maxIter)
    {
      NumericVector wOld = clone(w); 
      NumericVector wNew(k);
      
      NumericVector term1(k, 0.0); 
      double term2 = 0.0; 
      
      for (int i = 0; i < n; i++)
      {
        double wTz = 0.0; 
        for (int j = 0; j < k; j++) wTz += w[j] * Z(i, j); 
        
        double gVal = gFunc(wTz); 
        double gPrimeValue = gPrime(wTz);
        
        for (int j = 0; j < k; j++) term1[j] += Z(i, j) * gVal; 
        term2 += gPrimeValue;
      }
      
      // Gram-Schmidt
      for (int j = 0; j < c; j++)
      {
        double dot = 0.0; 
        for (int i = 0; i < k; i++) dot += wNew[i] * W(i, j); 
        for (int i = 0; i < k; i++) wNew[i] -= dot * W(i, j);
      }
      
      normalizeVector(wNew); 
      
      double convergenceCheck = 0.0; 
      for (int j = 0; j < k; j++) convergenceCheck += wNew[j] * wOld[j]; 
      
      w = wNew; 
      if (std::abs(std::abs(convergenceCheck) - 1.0) < tol)
      {
        converged = true; 
      }
      iter++; 
    }
    for (int j = 0; j < k; j++) W(j, c) = w[j]; 
  }
  // Projection
  NumericMatrix S = matrixMultiply(Z, W); 
  NumericMatrix WFull = matrixMultiply(Vk, W); 
  
  return List::create(
    Named("S") = S,                 
    Named("Wwhitened") = W,        
    Named("Wunmixing") = WFull,   
    Named("WhitenedData") = Z      
  );
}
