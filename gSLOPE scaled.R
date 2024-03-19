library("grpSLOPE")
library(glasso)

proximity_matrix_standardized <- function(matrix_in, lambdas) {
  output <- matrix_in
  precision_entries <- matrix_in[lower.tri(matrix_in, FALSE)]
  calculated_entries <- grpSLOPE::prox_sorted_L1(as.matrix(precision_entries),
    lambdas,
    method = c("c")
  )
  output[lower.tri(output, FALSE)] <- calculated_entries
  output[upper.tri(output, FALSE)] <- calculated_entries
  diag(output) <- 1.

  output
}


gslope_standardized_new <- function(
    sample_cov,
    lambdas,
    rho = 1.1,
    max_iter = 500,
    epsilon = 1e-04,
    progress = TRUE) {
  if (!(nrow(sample_cov) == ncol(sample_cov))) {
    stop("Covariance matrix must be square.")
  }

  # Parameters initialization
  Z <- sample_cov # Initialized to zero, probably it is the best choice
  Y <- Z
  diag(Y) <- 1
  X <- diag(nrow(sample_cov))

  # Start iteration
  if (progress) {
    for (iter in 1:max_iter) {
      C_tilde <- Y - Z - (sample_cov / rho)

      # Perform the eigenvalue decomposition
      C_eigen <- eigen(C_tilde, symmetric = TRUE)
      C_eigen_val <- C_eigen$val # Eigenvalues
      C_eigen_vec <- C_eigen$vec # Eigenvectors

      # Formula implementation
      F_rho <- 0.5 * diag(C_eigen_val + sqrt(C_eigen_val^2 + 4 / rho))
      X <- C_eigen_vec %*% F_rho %*% t(C_eigen_vec)

      Y_old <- Y
      Y <- proximity_matrix_standardized(X + Z, lambdas / rho) # CHANGE FROM proximity_matrix

      # Update step
      Z <- Z + rho * (X - Y)

      # Compute the primal and dual gap
      primal_residual <- norm(X - Y, type = "F")
      dual_residual <- norm(rho * (Y - Y_old), type = "F")

      #    print(paste0("Current progress Gslope: ", (iter / max_iter) * 100))

      # Stop condition
      if (primal_residual < epsilon & dual_residual < epsilon) {
        break
      }
    }
  } else {
    for (iter in 1:max_iter) {
      C_tilde <- Y - Z - (sample_cov / rho)

      # Perform the eigenvalue decomposition
      C_eigen <- eigen(C_tilde, symmetric = TRUE)
      C_eigen_val <- C_eigen$val # Eigenvalues
      C_eigen_vec <- C_eigen$vec # Eigenvectors

      # Formula implementation
      F_rho <- 0.5 * diag(C_eigen_val + sqrt(C_eigen_val^2 + 4 / rho))
      X <- C_eigen_vec %*% F_rho %*% t(C_eigen_vec)

      Y_old <- Y
      Y <- proximity_matrix_standardized(X + Z, lambdas / rho) # change from proximity_matrix

      # Update step
      Z <- Z + rho * (X - Y)

      # Compute the primal and dual gap
      primal_residual <- norm(X - Y, type = "F")
      dual_residual <- norm(rho * (Y - Y_old), type = "F")

      # Stop condition
      if (primal_residual < epsilon & dual_residual < epsilon) {
        break
      }
    }
  }

  X[abs(X) < 1e-04] <- 0 # Thresholding if abs(entries) <= 10e-4

  list(
    precision_matrix = X,
    cov_matrix = solve(X),
    iterations = iter,
    prim_res = primal_residual,
    dual_res = dual_residual
  )
}


###################### BK functions

myNorm <- function(x, a) {
  e <- rep(1, ncol(a))
  return(norm(diag(x) %*% a %*% diag(x) %*% e - e, type = "2"))
}

check <- function(a, max_iter = 500) {
  n <- ncol(a)
  e <- rep(1, n)
  if (myNorm(e, a) < 1) {
    return(as.vector(e))
  }
  if (myNorm(1 / diag(a), a) < 1) {
    return(as.vector(1 / diag(a)))
  }

  iter <- 0
  b <- e - a %*% e
  t0 <- 1
  x0 <- e

  x0 <- x0 + one_step(a, x0, b, t0)
  hatx <- sqrt(t0) * x0

  while (myNorm(hatx, a) > 3 / 4 && iter < max_iter) {
    t0 <- t0 * (1 - 1 / (4 * sqrt(n) + 1))
    x0 <- x0 + one_step(a, x0, b, t0)

    hatx <- sqrt(t0) * x0
    iter <- iter + 1
  }

  return(as.vector(hatx))
}
#########################################
one_step <- function(a, x0, b = 0, t = 1) { # a matrix, x0 vector
  rhs <- 1 / x0 - t * a %*% x0 - t * b
  lhs <- diag(x0^(-2)) + a
  y <- backsolve(lhs, rhs)
  return(as.vector(y))
}
#########################################
many_steps <- function(a, x0, tol = 1e-8, max_iter = 500) {
  iter <- 0
  x <- x0
  correction <- one_step(a, x)
  while (max(abs(correction)) > tol && iter < max_iter) {
    x <- x + correction
    iter <- iter + 1
    correction <- one_step(a, x)
  }
  return(x + correction)
}
#########################################
solve_for_D <- function(RS) {
  x0 <- check(RS)
  Dvec <- many_steps(RS, x0)
  D <- diag(Dvec)
  return(D)
}
#########################################
pcSLOPE <- function(lambda, S, tol = 1e-6, max_iter = 100) {
  diagS <- diag(S)
  S <- cov2cor(S)
  p <- ncol(S)
  iter <- 0
  INITIAL <- gslope_standardized_new(S, lambda)
  R <- INITIAL$precision_matrix
  Ri <- INITIAL$cov_matrix
  error <- 1
  D <- diag(p)

  while (error > tol && iter < max_iter) {
    Rold <- R
    Riold <- Ri
    D <- solve_for_D(Rold * S)
    inputMatrix <- D %*% S %*% D - D * D + diag(diag(Riold))
    # mineigen = min(eigen(inputMatrix)$values)
    #   print(mineigen)
    # if(mineigen<0 && -mineigen>=lambda){
    #  print("ERROR2")
    #  return(FALSE)
    # }
    GSLOPE <- gslope_standardized_new(inputMatrix, lambda)
    R <- GSLOPE$precision_matrix
    Ri <- GSLOPE$cov_matrix

    error <- max(abs(R - Rold))

    iter <- iter + 1
  }
  D <- D %*% diag(diagS^(-1 / 2))
  K <- D %*% R %*% D
  myList <- list(R, Ri, D, K, iter)

  #  print(iter)

  names(myList) <- c("R", "Ri", "D", "K", "iter")
  return(myList)
}

########### SIMULATIONS###########
########### SIMULATIONS###########
########### SIMULATIONS###########

p <- 10
n <- 100

val <- 1 / 2

trueK <- matrix(0, p, p)
for (i in 1:p) {
  trueK[i, i] <- 1
  if (i > 1) {
    trueK[i, i - 1] <- val
  }
  if (i < p) {
    trueK[i, i + 1] <- val
  }
}
trueSigma <- solve(trueK)
trueR <- trueK

Z <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = trueSigma)
S <- (t(Z) %*% Z) / n


lambda2 <- seq(p * (p - 1) / 2, 1) / (p * (p - 1) / 2)
pcSLOPE(lambda2 * 0.5, S)
