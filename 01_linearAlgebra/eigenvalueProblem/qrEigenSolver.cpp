#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// ===============================================================================
// Utility mathematical functions 
// ===============================================================================
// QR decomposition
List QRDecompose(NumericMatrix mat){
  int m = mat.nrow(); // rows
  int n = mat.ncol(); // cols
  
  NumericMatrix Q = clone(mat); 
  NumericMatrix R(n, n); 
  
  for(int k = 0; k < n; k++){
    double norm_sq = 0.0; 
    for(int i = 0; i < m; i++){
      norm_sq += Q(i, k)*Q(i, k); 
    }
    R(k, k) = std::sqrt(norm_sq); 
    
    for(int i = 0; i < m; i++){
      Q(i, k) /= R(k, k); 
    }
    
    for(int j = k + 1; j < n; j++){
      double dot = 0.0; 
      for(int i = 0; i < m; i++){
        dot += Q(i, k)*Q(i, j); 
      }
      R(k, j) = dot; 
      
      for(int i = 0; i < m; i++){
        Q(i, j) -= R(k, j) * Q(i, k);
      }
    }
  }
  
  return List::create(
    Named("Q") = Q, 
    Named("R") = R
  ); 
}

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

// build the identity matrix
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
    Named("eigenvectors") = Qtotal, 
    Named("Ak") = Ak
  );
}
