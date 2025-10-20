n_zero_eval = function(eval, th = 1e-10){
  # returns number of zero eigen values
  return(sum(eval < th))
}

eval_zero = function(eval, k){
  eval_res = eval;
  eval_res[1:k] = 0.0
  return(eval_res)
}

GS_complete = function(Evec, k){
  n = dim(Evec)[1]
  if (k < n){
    Q = matrix(NA, n, k)
    X = Evec[,1:k]
    X[,1] = rep(1/sqrt(n), n)
    Q[,1] = rep(1/sqrt(n), n)
    V = X 
    for (i in 2:k){
      for (j in 1:(i-1)){
        V[,i] = V[,i] -((t(Q[,j])%*%X[,i])%*%Q[,j])
      }
      Q[,i] = V[,i]/sqrt(sum(V[,i]*V[,i]))
    }
    return(cbind(Q, Evec[,(k+1):n]))
  }else {
    stop("the number of zero eigenvalue is", k, "out of", n)
  }
}

cen_eigen = function(K){
  n = dim(K)[1]
  E = eigen(K)
  l = E$values[n:1]
  Q = E$vectors[,n:1]
  k = n_zero_eval(l)
  if(k>0){
    l = eval_zero(l,k)
    if(k>1){
      Q = GS_complete(Q,k)
    } else {
      Q[,1] = rep(1/sqrt(n), n) 
    }
    return(list('values'= l,'vectors'= Q)) 
  } else{
    stop('k must be positive; found k=' , as.character(k))
  }
}

# mat_vec_prod = function(A,s, nd, nnd){
#   S = matrix(s, nd,nnd)
#   return(c(t(A%*%S)))
# }

mat_vec_prod = function (A,s){
  nd = ncol(A)
  nnd = length(s)/nd
  S = matrix(s, nd,nnd)
  return(c(t(A%*%S)))
}

mat_vec_prod_2d <- function(A, B, u) {
  ncolA <- ncol(A)
  ncolB <- ncol(B)
  Umat <- matrix(u, nrow = ncolB, ncol = ncolA)
  result <- B %*% Umat %*% t(A)
  as.vector(result)
}

kronecker_matvec <- function(A_list, s) {
  stopifnot(length(A_list) >= 2)
  res <- s
  for (i in seq(length(A_list), 2, by = -1)) {
    res <- mat_vec_prod_2d(A_list[[i-1]], A_list[[i]], res)
  }
  return(res)
}

mat_mat_prod <- function(A, B, C) {
  ncolA <- ncol(A)
  ncolB <- ncol(B)
  M <- matrix(0, nrow = nrow(A) * nrow(B), ncol = ncol(C))
  for (i in seq_len(ncol(C))) {
    Cmat <- matrix(C[, i], nrow = ncolB, ncol = ncolA)
    M[, i] <- as.vector(B %*% (Cmat %*% t(A)))
  }
  M
}

compute_post_var_update_3d <- function(A1,A2,A3,b){
  # computing (A1\otimesA2\otimesA3)'(diag(b))^{-1}(A1\otimesA2\otimesA3)'
  p1 <- ncol(A1)
  p2 <- ncol(A2)
  p3 <- ncol(A3)
  P <- p1 * p2 * p3
  N <- length(b)
  
  M <- matrix(0, nrow = P, ncol = P)
  
  for (u in 1:P) {
    eu <- rep(0, P)
    eu[u] <- 1
    ku <- mat_mat_prod(A1, diag(1, nrow(A2)), mat_mat_prod(A2, A3, matrix(eu, ncol=1)))
    wu <- ku / b
    
    for (v in u:P) {
      ev <- rep(0, P)
      ev[v] <- 1
      kv <- mat_mat_prod(A1, diag(1, nrow(A2)), mat_mat_prod(A2, A3, matrix(ev, ncol=1)))
      M[u, v] <- sum(wu * kv)
      if (u != v) M[v, u] <- M[u, v]
    }
  }
  return(M)
}

compute_post_var_update_3d_any_rank <- function(A1, A2, A3, lambda) {
  N1 <- nrow(A1)
  N2 <- nrow(A2)
  N3 <- nrow(A3)
  
  stopifnot(length(lambda) == N1 * N2 * N3)
  
  # Reshape lambda
  Larray <- array(lambda, dim = c(N1, N2, N3))
  
  # Compute marginal weights
  W1 <- apply(1 / Larray, 1, sum)  # sum over j,k
  W2 <- apply(1 / Larray, 2, sum)  # sum over i,k
  W3 <- apply(1 / Larray, 3, sum)  # sum over i,j
  
  # Handle each Ai
  get_component <- function(Ai, Wi) {
    if (ncol(Ai) == 1) {
      sum(Wi)
    } else {
      t(Ai) %*% (Wi * Ai)
    }
  }
  
  M1 <- get_component(A1, W1)
  M2 <- get_component(A2, W2)
  M3 <- get_component(A3, W3)
  
  # Final Kronecker
  if (length(M1) == 1 && length(M2) == 1 && length(M3) == 1) {
    M <- M1 * M2 * M3
  } else if (length(M1) == 1 && length(M2) == 1) {
    M <- M1 * M2 * M3
  } else if (length(M1) == 1 && length(M3) == 1) {
    M <- M1 * kronecker(M2, M3)
  } else if (length(M2) == 1 && length(M3) == 1) {
    M <- M2 * M3 * M1
  } else if (length(M1) == 1) {
    M <- kronecker(M2, M3) * M1
  } else if (length(M2) == 1) {
    M <- kronecker(M1, M3) * M2
  } else if (length(M3) == 1) {
    M <- kronecker(M1, M2) * M3
  } else {
    M <- kronecker(kronecker(M1, M2), M3)
  }
  
  return(M)
}

compute_post_var_diag_3d <- function(A1, A2, A3, lambda) {
  N1 <- nrow(A1)
  N2 <- nrow(A2)
  N3 <- nrow(A3)
  
  p1 <- ncol(A1)
  p2 <- ncol(A2)
  p3 <- ncol(A3)
  
  stopifnot(length(lambda) == N1 * N2 * N3)
  
  # Reshape lambda
  Larray <- array(lambda, dim = c(N1, N2, N3))
  Warray <- 1 / Larray
  
  diagM <- numeric(p1 * p2 * p3)
  idx <- 1
  
  for (i1 in 1:p1) {
    A1sq <- A1[, i1]^2
    for (i2 in 1:p2) {
      A2sq <- A2[, i2]^2
      for (i3 in 1:p3) {
        A3sq <- A3[, i3]^2
        
        # Element-wise multiply
        term <- outer(A1sq, outer(A2sq, A3sq, "*"), "*")
        diagM[idx] <- sum(term * Warray)
        idx <- idx + 1
      }
    }
  }
  return(diagM)
}

pcg = function(Q1,Q2,eval,d,y,initx, tol = 1e-3, itmax=1000){
  # output : x = (K + D)^{-1} y  with preconditioning matrix C = D^{-1/2} and D = diag_matrix(d)
  #          note K = Q %*% diag(eval) %*% Q' where Q = Q1 \otime Q2
  n1 = dim(Q1)[1] 
  n2 = dim(Q2)[1]
  n = n1*n2
  #dinv = (1/d) 
  #initial values
  x = initx
  v = eval*kron_mat_vec_2d(t(Q1), t(Q2), x) # diag_matrix(eval)%*%(Q)%*%x
  Kx = kron_mat_vec_2d(Q1,Q2,v)             # Q%*%v
  r = y - Kx - (d*x)                        # y- (K+D)x = y - Kx - Dx
  z = r/d                                   # D^{-1}%*%r , preconditioning
  p = z
  k = 0
  while (k<itmax){
    rz_old = c(crossprod(r,z))
    #  p'.(K+D).p = p'.(K.p + D.p)
    v = eval*kron_mat_vec_2d(t(Q1), t(Q2), p)
    KaddD.p = kron_mat_vec_2d(Q1,Q2,v) + d*p # (K+D)p = Kp + Dp
    alpha = rz_old/c(crossprod(p, KaddD.p))
    x = x + alpha*p
    r = r - alpha*(KaddD.p)
    if ( (sqrt(c(crossprod(r))) < tol) == TRUE){
      break
    } 
    z = r/d  #D^{-1}%*%r 
    beta = c(crossprod(r,z))/rz_old
    p = z + beta*p
    k = k + 1
  }
  return(list('sol'= x, 'iteration'=k))
}
