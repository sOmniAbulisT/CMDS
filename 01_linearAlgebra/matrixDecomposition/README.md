# Matrix Decomposition

This module implements the three core matrix decomposition algorithms in numerical linear algebra from scratch using C++ (`Rcpp`). 
The project aims to bypass high-level language encapsulation to deeply explore the underlying mathematical logic, memory manipulation, 
and numerical stability in floating-point operations.

## Core Implementations
This directory contains three independent C++ modules, each designed for specific numerical computing scenarios:

1. LU Decomposition (`luDecomposition.cpp`)

* Mathematical Principle: Factorizes a square matrix into the product of a lower triangular matrix and an upper 
  triangular matrix, denoted as $PA = LU$.  
  
* Implementation Highlights: Utilizes Doolittle's Method with a rigorous implementation of Partial Pivoting. A permutation matrix 
  $P$ is constructed to record row swaps, ensuring exceptional numerical stability even when pivot elements are close to zero.
  
* Applications: Solving systems of linear equations, matrix inversion, and determinant calculation.

2. Cholesky Decomposition (`choleskyDecomposition.cpp`)

* Mathematical Principle: Factorizes a Symmetric Positive-Definite matrix into $A = LL^T$.

* Implementation Highlights: Highly optimized for covariance matrices. During backward substitution and 
  inversion, direct index swapping is leveraged instead of explicit matrix transposition in memory, significantly reducing space 
  complexity.
  
* Applications: Multivariate normal sampling, kernel matrix operations in Gaussian Processes and Bayesian Optimization.

3. QR Decomposition (`qrDecomposition.cpp`)

* Mathematical Principle: Factorizes a matrix into an orthogonal matrix and an upper triangular matrix, denoted as $A = QR$.

* Implementation Highlights: Abandons Classical Gram-Schmidt in favor of the Modified Gram-Schmidt orthogonalization process. 
  This architecture drastically reduces the accumulation of round-off errors in floating-point arithmetic, guaranteeing strict 
  orthogonality of the $Q$ matrix ($Q^TQ = I$).
  
* Applications: Solving overdetermined systems, highly stable Ordinary Least Squares regression solvers.