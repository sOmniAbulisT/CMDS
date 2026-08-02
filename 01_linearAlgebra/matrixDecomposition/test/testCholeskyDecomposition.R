library(Rcpp)
library(microbenchmark)

# 載入 C++ 程式碼
sourceCpp("01_linearAlgebraicTechniques/matrixDecomposition/src/choleskyDecomposition.cpp")

# 1. 建立測試資料 (必須是對稱正定矩陣)
set.seed(42)
n <- 100
# 產生隨機矩陣 Z，並透過 Z^T * Z 構造對稱正定矩陣 A
Z <- matrix(rnorm(n * n), nrow = n, ncol = n)
A <- crossprod(Z) + diag(1e-5, n) # 加上微小對角線元素確保嚴格正定
b <- rnorm(n)

# 2. 測試 Cholesky 分解正確性
cat("=== 測試 CholDecompose ===\n")
chol_res <- CholDecompose(A)
L <- chol_res$L
# 驗證 A = L * L^T
LLT <- L %*% t(L)
error_Chol <- max(abs(A - LLT))
cat("A 與 L*L^T 的最大絕對誤差: ", error_Chol, "\n\n")

# 3. 測試線性方程求解 (Ax = b)
cat("=== 測試 CholSolve ===\n")
x_custom <- CholSolve(chol_res, b)
x_base <- solve(A, b) 
error_solve <- max(abs(x_custom - x_base))
cat("求解 Ax=b 的最大絕對誤差: ", error_solve, "\n\n")

# 4. 測試反矩陣
cat("=== 測試 CholInvert ===\n")
inv_custom <- CholInvert(chol_res)
inv_base <- solve(A)
error_invert <- max(abs(inv_custom - inv_base))
cat("反矩陣的最大絕對誤差: ", error_invert, "\n\n")

# 5. 測試行列式
cat("=== 測試 CholDeterminant ===\n")
det_custom <- CholDeterminant(chol_res)
det_base <- det(A)
cat("自訂行列式: ", det_custom, "\n")
cat("R 內建行列式: ", det_base, "\n")
cat("相對誤差: ", abs(det_custom - det_base) / abs(det_base), "\n\n")