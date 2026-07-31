library(Rcpp)
library(microbenchmark)

sourceCpp("01_linearAlgebraicTechniques/matrixDecomposition/src/luDecomposition.cpp")

# 1. test data
set.seed(42)
n <- 100

A <- matrix(rnorm(n * n), nrow = n, ncol = n)
b <- rnorm(n)

# 2. 測試 LU 分解正確性
cat("=== 測試 LU 分解 ===\n")
lu_res <- LUDecompose(A)
# 驗證 PA = LU
PA <- lu_res$P %*% A
LU <- lu_res$L %*% lu_res$U
# 計算誤差，應該要極度接近 0
error_LU <- max(abs(PA - LU))
cat("PA 與 LU 的最大絕對誤差: ", error_LU, "\n\n")

# 3. 測試線性方程求解 (Ax = b)
cat("=== 測試 LUSolve ===\n")
x_custom <- LUSolve(lu_res, b)
x_base <- solve(A, b) # R 內建的求解函數
error_solve <- max(abs(x_custom - x_base))
cat("求解 Ax=b 的最大絕對誤差: ", error_solve, "\n\n")

# 4. 測試反矩陣
cat("=== 測試 LUInvert ===\n")
inv_custom <- LUInvert(lu_res)
inv_base <- solve(A)
error_invert <- max(abs(inv_custom - inv_base))
cat("反矩陣的最大絕對誤差: ", error_invert, "\n\n")

# 5. 測試行列式
cat("=== 測試 LUDeterminant ===\n")
det_custom <- LUDeterminant(lu_res)
det_base <- det(A)
cat("自訂行列式: ", det_custom, "\n")
cat("R 內建行列式: ", det_base, "\n")
cat("相對誤差: ", abs(det_custom - det_base) / abs(det_base), "\n\n")

# 6. 簡單效能評估 (Benchmark)
cat("=== 求解效能比較 ===\n")
benchmark_res <- microbenchmark(
  Custom_LU_Solve = LUSolve(LUDecompose(A), b),
  Base_R_Solve = solve(A, b),
  times = 100
)
print(benchmark_res)
