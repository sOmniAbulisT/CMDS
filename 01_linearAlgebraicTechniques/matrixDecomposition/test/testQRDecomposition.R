library(Rcpp)
library(microbenchmark)

# 載入 C++ 程式碼
sourceCpp("01_linearAlgebraicTechniques/matrixDecomposition/src/qrDecomposition.cpp")

# 1. 建立測試資料
set.seed(42)
m <- 100 # 列數 (觀測值)
n <- 50  # 行數 (特徵數) - 這裡刻意設計 m > n 來測試非方陣

# 產生隨機矩陣 A (m x n) 與向量 b (m x 1)
A <- matrix(rnorm(m * n), nrow = m, ncol = n)
b <- rnorm(m)

# 2. 測試 QR 分解正確性
cat("=== 測試 QRDecompose (非方陣) ===\n")
qr_res <- QRDecompose(A)
Q <- qr_res$Q
R <- qr_res$R

# 驗證 A = QR
QR_reconstruct <- Q %*% R
error_QR <- max(abs(A - QR_reconstruct))
cat("A 與 Q*R 的最大絕對誤差: ", error_QR, "\n")

# 驗證 Q 的正交性 (Q^T * Q = I)
QTQ <- t(Q) %*% Q
error_ortho <- max(abs(QTQ - diag(1, n)))
cat("Q^T * Q 與單位矩陣 I 的最大絕對誤差: ", error_ortho, "\n\n")

# 3. 測試線性方程求解 (最小平方法 OLS)
cat("=== 測試 QRSolve (最小平方法) ===\n")
x_custom <- QRSolve(qr_res, b)
# R 內建的最小平方法求解 (lm.fit 或 qr.solve)
x_base <- qr.solve(A, b) 
error_solve <- max(abs(x_custom - x_base))
cat("求解 OLS 的最大絕對誤差: ", error_solve, "\n\n")

# ==========================================
# 針對方陣的測試 (Inverse 與 Determinant)
# ==========================================
cat("=== 方陣專屬測試 (Inverse & Determinant) ===\n")
A_sq <- matrix(rnorm(n * n), nrow = n, ncol = n)
qr_sq_res <- QRDecompose(A_sq)

# 測試反矩陣
inv_custom <- QRInvert(qr_sq_res)
inv_base <- solve(A_sq)
error_invert <- max(abs(inv_custom - inv_base))
cat("反矩陣的最大絕對誤差: ", error_invert, "\n")

# 測試行列式 (驗證絕對值)
det_custom <- QRDeterminant(qr_sq_res)
det_base <- det(A_sq)
cat("自訂行列式: ", det_custom, "\n")
cat("R 內建行列式 (取絕對值): ", abs(det_base), "\n")
cat("相對誤差: ", abs(det_custom - abs(det_base)) / abs(det_base), "\n\n")