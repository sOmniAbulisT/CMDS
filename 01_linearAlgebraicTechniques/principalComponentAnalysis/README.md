#Principal Component Analysis
A high-performance, dependency-free implementation of Principal Component Analysis (PCA) built from scratch using C++ and 
seamlessly integrated into R via `Rcpp`.

This directory demonstrates the fundamental numerical linear algebra behind dimensionality reduction. Instead of relying on standard 
black-box linear algebra libraries (like LAPACK or BLAS), this engine implements a custom QR Algorithm to solve the Eigenvalue problem, 
showcasing a deep understanding of algorithmic stability and matrix transformations.

##Core Features
* Zero External Math Dependencies: Eigenvalues and eigenvectors are computed using a custom-built QR Algorithm with 
  Modified Gram-Schmidt orthogonalization.
  
* Memory-Aware Centering: Implements column-major data centering to maximize CPU cache hit rates during large-scale matrix operations.

* Automated Feature Selection: Dynamically sorts eigenvalues and their corresponding eigenvectors in descending order, ensuring accurate principal 
  component selection.
  
* Rcpp Integration: Combines the extreme computational speed of C++ with the statistical ecosystem of R.

##Directory Structure
```
01_linearAlgebraicTechniques/
├── src/
│   └── pca.cpp            # The main PCA pipeline (Centering, Covariance, Projection)
├── testPCA.R             # Unit tests comparing custom PCA against R's native prcomp()
└── README.md
```

##The 5-Step Algorithm Pepline
The `PCA` function in `src/pca.cpp` strictly follows the standard mathematical approach for Principal Component Analysis:

1. Data Centering: Computes the mean of each feature (column) and subtracts it from the respective observations.

2. Covariance Matrix Computation: Calculates the covariance matrix $\Sigma = \frac{1}{n-1} X_c^T X_c$, utilizing matrix symmetry to 
  optimize calculation time.
  
3. Solve Eigenvalue Problem: Calls the custom `QREigenSolver` to diagonalize the covariance matrix via similarity transformations 
  ($A_k = R_k Q_k$).
  
4. Feature Vector Selection: Sorts the resulting eigenvalues in descending order and dynamically reorders the corresponding eigenvector columns 
  (Loadings).
  
5. Data Recast (Projection): Projects the centered original data onto the selected principal component axes to produce the final dimension-reduced coordinates 
  (PC Scores).
  
##Getting Started

###Prerequisites

* R environment
* C++ compiler (GCC/Clang)
* `Rcpp` package installed in R

###Execution & Testing
Run the provided R script to compile the C++ source code dynamically and verify the precision against R's native `prcomp()` function.

```R
# Install Rcpp if you haven't already
# install.packages("Rcpp")

library(Rcpp)

# Compile and load the C++ PCA engine
sourceCpp("01_linearAlgebraicTechniques/principalComponentAnalysis/src/pca.cpp")

# Generate a highly-correlated synthetic dataset
set.seed(42)
n <- 100
p <- 5
X <- matrix(rnorm(n * p), nrow = n)
X[, 2] <- X[, 1] * 2.5 + rnorm(n, sd = 0.5)

# Run the custom PCA, extracting the top 2 Principal Components
my_pca <- RunPCA(X, k = 2)

# View the extracted Principal Component Scores
print(head(my_pca$PC_scores))
```