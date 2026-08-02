#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

/**
 * QR decomposition using Modified Gram-Schmidt
 */
// [[Rcpp::export]]
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

/**
 * Solve the system of linear equation (Ax=b) using QR decomposition
 * 
 * 
 */

//[[Rcpp::export]]
NumericVector QRSolve(List qr_result, NumericVector b){
  NumericMatrix Q = qr_result["Q"]; 
  NumericMatrix R = qr_result["R"]; 
  
  int m = Q.nrow(); // row
  int n = Q.ncol(); // col
  if(b.size() != m) stop("Vector b dimension must match the number of rows in Q. "); 
  
  NumericVector y(n); 
  for(int i = 0; i < n; i++){
    double sum = 0.0; 
    for(int j = 0; j < m; j++){
      sum += Q(j, i) * b[j];
    }
    y[i] = sum; 
  }
  
  NumericVector x(n); 
  for(int i = n-1; i >= 0; i--){
    double sum = 0.0; 
    for(int j = i+1; j < n; j++){
      sum += R(i, j)*x[j]; 
    }
    x[i] = (y[i] - sum) / R(i, i); 
  }
  
  return x; 
}

/**
 * Calculate the Inverse Matrix using QR Decomposition
 * 
 */
// [[Rcpp::export]]
NumericMatrix QRInvert(List qr_result){
  NumericMatrix Q = qr_result["Q"]; 
  NumericMatrix R = qr_result["R"];
  
  int m = Q.nrow(); //row
  int n = Q.ncol();
  if(m != n) stop("Matrix must be square. ");
  
  NumericMatrix IA(n, n); 
  
  for(int j = 0; j < n; j++){
    NumericVector y(n);
    for(int i = 0; i < n; i++){
      y[i] = Q(j, i);
    }
    
    for(int i = n - 1; i >= 0; i--){
      double sum = 0.0; 
      for(int k = i + 1; k < n; k++){
        sum += R(i, k)*IA(k, j); 
      }
      IA(i, j) = (y[i] - sum) / R(i, i); 
    }
  }
  
  return IA; 
}

/**
 * Calculate the Determinant of a Matrix using QR Decomposition
 */

//[[Rcpp::export]]
double QRDeterminant(List qr_result){
  NumericMatrix R = qr_result["R"]; 
  
  int m = R.nrow(); 
  int n = R.ncol(); 
  if(m != n) stop("Matrix must be square"); 
  
  double det = 1.0;
  for(int i = 0; i < n; i++){
    det *= R(i, i); 
  }
  
  return det;
}