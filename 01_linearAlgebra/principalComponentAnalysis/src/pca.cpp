#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace Rcpp;

#include "01_linearAlgebraicTechniques/eigenvalueProblem/qrEigenSolver.cpp"

/**
 * Principal Component Analysis
 * @param X: input data matrix
 * @param k: dimension reduction target 
 */
//[[Rcpp::export]]
List PCA (NumericMatrix X, int k)
{
  int n = X.nrow(); 
  int p = X.ncol(); 
  
  if (k > p) stop("k must be less than or equal to p"); 
  
  // ====================
  // Step 1. Center data
  // ====================
  NumericMatrix Xc = clone(X);
  for (int j = 0; j < p; j++)
  {
    double sum = 0.0; 
    for (int i = 0; i < n; i++) sum += Xc(i, j); 
    double mean = sum / n; 
    
    for (int i = 0; i < n; i++) Xc(i, j) -= mean; 
  }
  
  // ======================================
  // Step 2. Compute the covariance matrix
  // ======================================
  NumericMatrix Cov(p, p); 
  for (int i = 0; i < p; i++)
  {
    for (int j = i; j < p; j++) // symmeric compute upper triangular
    {
      double dotProduct = 0.0; 
      for (int r = 0; r < n; r++)
      {
        dotProduct += Xc(r, i) * Xc(r, j);
      }
      Cov(i, j) = dotProduct / (n - 1);
      Cov(j, i) = Cov(i, j); // copy to lower triangular
    }
  }
  
  // =================================
  // Step 3. Solve eigenvalue problem
  // =================================
  List eigenRes = QREigenSolver(Cov); 
  NumericVector eigenvalues = eigenRes["eigenvalues"]; 
  NumericVector eigenvector = eigenRes["eigenvectors"]; 
  
  // ==================================
  // Step 4. Feature vector selection
  // ==================================
  // 1. bulid the index array
  std::vector<int> idx(p); 
  std::iota(idx.begin(), idx.end(), 0); 
  
  // 2. sort
  std::sort(idx.begin(), idx.end(), [&](int i, int j)
  {
    return eigenvalues[i] > eigenvalues[j]; 
  }); 
  
  // 3. 
  NumericMatrix vFeature(p, k); 
  NumericMatrix topEigenvalues(k); 
  
  for (int j = 0; j < k; j++)
  {
    int originalIdx = idx[j]; 
    topEigenvalues[j] = eigenvalues[originalIdx]; 
    for (int i = 0; i < p; i++)
    {
      vFeature(i, j) = eigenvectors(i, originalIdx);  
    }
  }
  
  // ============================
  // Step 5: Recast the data
  // ============================
  NumericMatrix Xfinal = matrixMultiply(Xc, vFeature); 
  
  return List::create(
    Named("PCscores") = Xfinal,       // (n x k)
    Named("loadings") = vFeature,      // (p x k)
    Named("eigenvalues") = topEigenvalues 
  );
}

