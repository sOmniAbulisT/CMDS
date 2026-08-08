#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;
#include "01_linearAlgebra/matrixDecomposition/src/qrDecomposition.cpp"

// matrix multiplication
NumericMatrix matrixMultiply (NumericMatrix A, NumericMatrix B)
{
  int n = A.nrow(); 
  NumericMatrix C(n, n); 
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      double sum = 0.0;
      for (int k = 0; k < n; k++)
      {
        sum += A(i, k) * B(k, j); 
      }
      C(i, j) = sum; 
    }
  }
  return C; 
}

// bulid the identity matrix
NumericMatrix identityMatrix (int n)
{
  NumericMatrix I(n, n); 
  for (int i = 0; i < n; i++) I(i, i) = 1.0; 
  return I; 
}

/**
 * QR Algorithm for Eigenvalues and Eigenvectors
 **/
//[[Rcpp::export]]
List QREigenSolver (NumericMatrix A, int maxIter = 200, double tol = 1e-10)
{
  int n = A.nrow();
  if (n != A.ncol()) stop("Matrix must be square");
  
  NumericMatrix Ak = clone(A); 
  NumericMatrix Qtotal = identityMatrix(n); 
  
  for (int k = 0; k < maxIter; k++)
  {
    List qrRes = QRDecompose(Ak); 
    NumericMatrix Q = qrRes["Q"]; 
    NumericMatrix R = qrRes["R"]; 
    
    NumericMatrix AkNew = matrixMultiply(R, Q); 
    Qtotal = matrixMultiply(Qtotal, Q); 
    
    double offDiagSum = 0.0; 
    for (int i = 1; i < n; i++)
    {
      for (int j = 0; j < i; j++)
      {
        offDiagSum += std::abs(AkNew(i, j)); 
      }
    }
    
    Ak = AkNew; 
    if (offDiagSum < tol) break; 
  }
  
  NumericVector eigenvalues(n);
  for (int i = 0; i < n; i++)
  {
    eigenvalues[i] = Ak(i, i);  
  }
  
  return List::create(
    Named("eigenvalues") = eigenvalues, 
    Named("eigenvectors") = Qtotal
  );
}
