#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace Rcpp;

#include "01_linearAlgebraicTechniques/eigenvalueProblem/qrEigenSolver.cpp"

// ==============================================================================
// Mathematical utilities
// ==============================================================================
// Transpose
NumericMatrix transposeMatrix (const NumericMatrix& A)
{
  int m = A.nrow(), n = A.ncol();
  NumericMatrix At(n, m); 
  
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      At(j, i) = A(i, j); 
    }
  }
  return At; 
}

// Multiplied
NumericMatrix multipliedMatrix (const NumericMatrix& A, const NumericMatrix& B)
{
  int m = A.nrow(), p = A.ncol(), n = B.ncol();
  NumericMatrix C(m, n); 
  
  for (int i = 0; i < m, i++)
  {
    for (int j = 0; j < n, j++)
    {
      double sum = 0.0; 
      for (int k = 0; k < p; k++)
      {
        sum += A(i, k) * B(k, j); 
      }
      C(i, j) = sum; 
    }
  }
  return C; 
}

// Normalization
void normalizeVector (NumericVector& v)
{
  double norm = 0.0; 
  for (int i = 0; i < v.size(); i++) norm += v[i] * v[i];
  norm = std::sqrt(norm);
  if (norm > 1e-12) 
  {
    for (int i = 0; i < v.size(); i++) v[i] /= norm;
  }
}

// ==============================================================================
// Singular Value Decomposition
// ==============================================================================
List SVD (NumericMatrix A, double tol = 1e-10)
{
  int n = A.nrow(); // number of row
  int m = A.ncol(); // number of column
  
  // ----------------------------------------------
  // Step 1: Compute AtA and eigen decomposition
  // ----------------------------------------------
  NumericMatrix At = transposeMatrix(A); 
  NumericMatrix AtA = multipliedMatrix(At, A);
  
  List eigenRes = QREigenSolver(AtA);
  NumericVector eigenvaluesRaw = eigenRes["eigenvalues"];
  NumericMatrix eigenvectorsRaw = eigenRes["eigenvectors"];
  
  // -----------------------------------------------
  // Step 2: Sort the eigen values
  // -----------------------------------------------
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](int i, int j) 
  {
    return eigenvaluesRaw[i] > eigenvaluesRaw[j];
  });
  
  // V matrix (n x n)
  NumericMatrix V(n, n);
  NumericVector singularValues(n);
  
  for (int j = 0; j < n; j++) 
  {
    int sortedIdx = idx[j];
    double lam = eigenvaluesRaw[sortedIdx];
    singularValues[j] = (lam > 0.0) ? std::sqrt(lam) : 0.0; 
    
    for (int i = 0; i < n; i++) 
    {
      V(i, j) = eigenvectorsRaw(i, sortedIdx);
    }
  }
  
  // S matrix (m x n diagonal matrix)
  NumericMatrix S(m, n); 
  int minDim = std::min(m, n); 
  for (int i = 0; i < minDim; i++)
  {
    S(i, i) = singularValues[i];
  }
  
  // ----------------------------------------------------
  // Step 3: Compute U matrix (m x m orthogonal matrix)
  // ----------------------------------------------------
  NumericMatrix U(m, m); 
  
  int rank = 0; 
  for (int j = 0; j < minDim; j++)
  {
    if (singularValues[j] > tol)
    {
      for (int i = 0; i < m; i++)
      {
        double sum = 0.0;
        for (int k = 0; k < n; k++) sum += A(i, k) * V(k, j);
        U(i, j) = sum / singularValues[j];
      }
      rank++; 
    } 
    else
    {
      break; 
    }
  }
  
  for (int j = rank; j < m; j++)
  {
    NumericVector uCand(m); 
    uCand[j % m] = 1.0; 
    
    for (int k = 0; k < j; k++)
    {
      double dot = 0.0; 
      for (int i = 0; i < m; i++) dot += uCand[i] * U(i, k); 
      for (int i = 0; i < m; i++) uCand[i] -= dot * U(i, k); 
    }
    
    normalizedVector(uCand); 
    
    double normCheck = 0.0;
    for (int i = 0; i < m; i++) normCheck += uCand[i] * uCand[i];
    if (normCheck < 1e-6) 
    {
      for (int i = 0; i < m; i++) uCand[i] = R::rnorm(0.0, 1.0);
      for (int k = 0; k < k < j; k++) 
      {
        double dot = 0.0;
        for (int i = 0; i < m; i++) dot += uCand[i] * U(i, k);
        for (int i = 0; i < m; i++) uCand[i] -= dot * U(i, k);
      }
      normalizeVector(uCand);
    }
    
    for (int i = 0; i < m; i++) U(i, j) = uCand[i];
  }
  
  return List::create(
    Named("U") = U,                    
    Named("S") = S,                     
    Named("V") = V,                     
    Named("d") = singularValues         
  );
}