## ---------- helpers ----------
squared_distance <- function(x, y) {
  dx <- x - y
  sum(dx * dx)
}
distance <- function(x, y) sqrt(squared_distance(x, y))
dot_self <- function(x) sum(x * x)

## ---------- kernels ----------
kernel_se <- function(x, y, rho) {
  exp(-squared_distance(x, y) / (2 * rho^2))
}
kernel_exponential <- function(x, y, rho) {
  exp(-distance(x, y) / rho)
}
kernel_matern32 <- function(x, y, rho) {
  r <- sqrt(3) * distance(x, y) / rho
  (1 + r) * exp(-r)
}
kernel_matern52 <- function(x, y, rho) {
  r <- sqrt(5) * distance(x, y) / rho
  (1 + r + (r^2) / 3) * exp(-r)
}
# fractional Brownian motion kernel (uses Hurst exponent H as 'param')
kernel_fBM <- function(x, y, H) {
  xysq <- squared_distance(x, y)
  xn2  <- dot_self(x)
  yn2  <- dot_self(y)
  # guard tiny negative due to FP error (pow of negative is undefined for non-integers)
  xysq <- max(xysq, 0)
  0.5 * (xn2^H + yn2^H - xysq^H)
}

## ---------- dispatcher (id as string: "se", "exp", "matern32", "matern52", "fbm") ----------
kernel_dispatch <- function(x, y, kernel, param) {
  kernel <- tolower(kernel)
  if (kernel == "se")          return(kernel_se(x, y, param))
  else if (kernel == "exp")    return(kernel_exponential(x, y, param))
  else if (kernel == "matern32") return(kernel_matern32(x, y, param))
  else if (kernel == "matern52") return(kernel_matern52(x, y, param))
  else if (kernel == "fbm")    return(kernel_fBM(x, y, param))
  stop("Unknown kernel. Use one of: 'se', 'exp', 'matern32', 'matern52', 'fbm'.")
}

## ----------  kernel matrix ----------
# - If Z is NULL (or identical to X): symmetric Gram with jitter on diag
# - Else: rectangular cross-kernel matrix K_ij = k(X_i, Z_j)
kernel_matrix <- function(X, Z = NULL, kernel, param, jitter = 0) {
  if (!is.matrix(X)) X <- as.matrix(X)
  
  # Symmetric Gram case
  if (is.null(Z) || identical(Z, X)) {
    n <- nrow(X)
    K <- matrix(0.0, n, n)
    for (i in seq_len(n)) {
      xi <- X[i, ]
      K[i, i] <- kernel_dispatch(xi, xi, kernel, param) + jitter
      if (i < n) {
        for (j in (i + 1L):n) {
          xj <- X[j, ]
          kij <- kernel_dispatch(xi, xj, kernel, param)
          K[i, j] <- kij
          K[j, i] <- kij
        }
      }
    }
    return(K)
  }
  
  # Rectangular cross-kernel case
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  if (ncol(X) != ncol(Z)) stop("X and Z must have the same number of columns.")
  n <- nrow(X); m <- nrow(Z)
  K <- matrix(NA_real_, n, m)
  for (i in seq_len(n)) {
    xi <- X[i, ]
    for (j in seq_len(m)) {
      zj <- Z[j, ]
      K[i, j] <- kernel_dispatch(xi, zj, kernel, param)
    }
  }
  return(K)
}

## ---------- centering a Gram matrix ----------
center_gram <- function(K) {
  n <- nrow(K)
  I <- diag(n)
  C <- I - matrix(1 / n, n, n)
  C %*% K %*% C
}
## ----------  centrerd kernel function ----------
# - kc(x,y) = k(x,y) - mean(k(x,X)) - mean(k(y,X)) + (1/n^2)sum(k(X,X))
# - not so efficient, better to use center_cross or center_self if possible
centred_kernel <- function(x, y = NULL, X, kernel, param) {
  n <- nrow(X)
  t1 <- kernel_dispatch(x, y, kernel, param)
  t2 <- kernel_matrix(x, X, kernel, param)
  if (is.null(y)) t3 <- t2 else t3 = kernel_matrix(y, X, kernel, param)
  t4 = kernel_matrix(X, kernel, param)
  return(t1 - mean(t2) - mean(t3) + (1/n^2)*sum(t4))
}

center_cross_mat <- function(Kcross_raw, Ktrain_raw) {
  N <- nrow(Ktrain_raw); n <- ncol(Kcross_raw)
  # column-wise means over training for each test point (length n), repeated down N rows
  term_cols <- matrix(colMeans(Kcross_raw), nrow = N, ncol = n, byrow = TRUE)
  # row means of training kernel, repeated across n columns
  term_rows <- matrix(rowMeans(Ktrain_raw), nrow = N, ncol = n)
  grand_mean_train <- sum(Ktrain_raw) / (N^2)
  Kcross_raw - term_cols - term_rows + matrix(grand_mean_train, N, n)
}

center_new_mat <- function(Knew_raw, Kcross_raw, Ktrain_raw) {
  n <- nrow(Knew_raw); N <- nrow(Ktrain_raw)
  # (1/N) 1_N^T Kcross_raw  => col means of cross: length n
  mZ <- colMeans(Kcross_raw)
  grand_mean_train <- sum(Ktrain_raw) / (N^2)
  Kc <- Knew_raw - matrix(mZ, n, n, byrow = TRUE) - matrix(mZ, n, n, byrow = FALSE) +
    matrix(grand_mean_train, n, n)
  # enforce symmetry numerically
  0.5 * (Kc + t(Kc))
}

## ---------- center cross-covector (new -> training) ----------
# kc(x*, X) = k(x*, X) - mean(k(x*, X))*1 - row_mean(K) + grand_mean(K)*1
center_cross <- function(k_raw, row_mean_train, grand_mean_train) {
  N <- length(k_raw)
  mean_new <- mean(k_raw)
  k_raw - mean_new * rep(1.0, N) - row_mean_train + grand_mean_train * rep(1.0, N)
}

## ---------- center self term ----------
# kc(x*, x*) = k(x*,x*) - 2*mean(k(x*,X)) + grand_mean(K)
center_self <- function(k_self_raw, k_raw, grand_mean_train) {
  mean_new <- mean(k_raw)
  k_self_raw - 2.0 * mean_new + grand_mean_train
}