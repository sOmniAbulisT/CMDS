#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp; 
using namespace arma; 

int findEnteringColumn (const mat& tableau)
{
  int lastRow = tableau.n_rows - 1; 
  int cols = tableau.n_cols - 1; 
  
  double minVal = 0.0; 
  int pivotCol = -1; 
  
  for (int j = 0; j < cols; ++j)
  {
    if (tableau(lastRow, j) < minVal)
    {
      minVal = tableau(lastRow, j); 
      pivotCol = j; 
    }
  }
  
  return pivotCol;
  
}

int performRatioTest (const mat& tableau, int pivotCol)
{
  int rows = tableau.n_rows - 1; 
  int lastCol = tableau.n_cols - 1; 
  
  double minRatio = datum::inf; 
  int pivotRow = -1; 
  
  for (int i = 0; i < rows; ++i)
  {
    double element = tableau(i, pivotCol); 
    
    if (element > 0)
    {
      double ratio = tableau(i, lastCol) / element; 
      if (ratio < minRatio)
      {
        minRatio = ratio; 
        pivotRow = i; 
      }
    }
  }
  
  return pivotRow; 
}

void executePivotStep (mat& tableau, int pivotRow, int pivotCol)
{
  double pivotVal = tableau(pivotRow, pivotCol);
  
  tableau.row(pivotRow) = tableau.row(pivotRow) / pivotVal; 
  
  for (uword i = 0; i < tableau.n_rows; ++i)
  {
    if (i != (uword)pivotRow)
    {
      double multiplier = tableau(i, pivotCol); 
      tableau.row(i) = tableau.row(i) - (multiplier * tableau.row(pivotRow)); 
    }
  }
}

//[[Rcpp::export]]
List runSimplexSolver(mat tableau)
{
  int maxIter = 1000; 
  int iter = 0; 
  
  while (iter < maxIter)
  {
    // Step 1
    int pivotCol = findEnteringColumn(tableau); 
    if (pivotCol == -1)
    {
      return List::create(Named("status") = "Optimal", 
                          Named("tableau") = tableau, 
                          Named("iterations") = iter); 
    }
    
    // Step 2
    int pivotRow = performRatioTest(tableau, pivotCol); 
    if (pivotRow == -1)
    {
      return List::create(Named("status") = "Optimal", 
                          Named("tableau") = tableau, 
                          Named("iterations") = iter); 
    }
    
    // Step 3
    executePivotStep(tableau, pivotRow, pivotCol); 
    iter++; 
  }
  
  return List::create(Named("status") = "Max Iterations Reached",
                      Named("tableau") = tableau);
}