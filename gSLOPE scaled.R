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

  while (myNorm(hatx, a) > 1 && iter < max_iter) {
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
  lhs <- diag(x0^(-2)) + t * a
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
# The funciton that coordinate descent is optimizing:
pcSLOPE_target_2x2 <- function(R, D, S, lambda) {
  my_trace <- function(M) {
    sum(diag(M))
  }
  r <- R[1, 2]
  log(1 - r * r) + 2 * log(det(D)) - my_trace(R %*% D %*% S %*% D) - 2 * lambda * abs(r)
}
#########################################
pcSLOPE_2x2 <- function(lambda, s, r, tol = 1e-6, max_iter = 100) {
  diagS <- rep(1, 2)
  S <- matrix(c(1, s, s, 1), nrow = 2)
  R <- matrix(c(1, r, r, 1), nrow = 2)

  pcSLOPE(lambda, S, R, tol, max_iter)
}
pcSLOPE <- function(lambda, S, R = NULL, tol = 1e-6, max_iter = 100) {
  diagS <- diag(S)
  S <- cov2cor(S)
  p <- ncol(S)
  iter <- 0
  Ri <- NULL
  if (is.null(R)) {
    INITIAL <- gslope_standardized_new(S, lambda)
    R <- INITIAL$precision_matrix
    Ri <- INITIAL$cov_matrix
  }
  error <- 1
  D <- diag(p)

  while (error > tol && iter < max_iter) {
    Rold <- R
    Riold <- Ri
    D_new <- solve_for_D(Rold * S)
    if (pcSLOPE_target_2x2(Rold, D_new, S, lambda) < -100) {
      stop("The matrix overflow. It used to make a singular D matrix. Investigate!")
    }
    D <- D_new
    inputMatrix <- D %*% S %*% D
    GSLOPE <- gslope_standardized_new(inputMatrix, lambda)
    R <- GSLOPE$precision_matrix
    diag(R) <- rep(1, p) # TODO: hacking; this should be true for output of gslope_standardized_new()
    Ri <- GSLOPE$cov_matrix

    error <- max(abs(R - Rold))

    iter <- iter + 1
  }
  D <- D %*% diag(diagS^(-1 / 2))
  K <- D %*% R %*% D
  myList <- list(R, Ri, D, K, iter, error)

  names(myList) <- c("R", "Ri", "D", "K", "iter", "error")
  return(myList)
}

########### SIMULATIONS###########

p <- 2
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


N_l <- 20
l_seq <- seq(0, 3 / 2, length.out = N_l + 2)[2:(N_l + 1)]
N_s <- 20
s_seq <- seq(-1, 1, length.out = N_s + 2)[2:(N_s + 1)]

args_grid <- tidyr::crossing(lambda = l_seq, s = s_seq)

single_r_result <- function(r, lambda, s) {
  pcSLOPE_res <- pcSLOPE_2x2(lambda, s, r)
  list(r_res = pcSLOPE_res$R[1, 2], n_iter = pcSLOPE_res$iter, error = pcSLOPE_res$error)
}

r_gap <- function(lambda, s) {
  N_r <- 20
  r_seq <- seq(-1, 1, length.out = N_r + 2)[2:(N_r + 1)]

  multiple_r_result <- sapply(r_seq, function(r) {
    single_r_result(r, lambda, s)$r_res
  })

  (max(multiple_r_result) - min(multiple_r_result))
}

start_time <- Sys.time()
result <- purrr::pmap_dbl(args_grid, r_gap) # 0.5 minut
max(result)
end_time <- Sys.time()
end_time - start_time

args_grid[which(result > 0.5), ]

lambda <- 1
s <- -0.905
N_r <- 1000

r_seq <- seq(-0.1, 0.9, length.out = N_r + 2)[2:(N_r + 1)]

multiple_r_result <- sapply(r_seq, function(r) {
  single_r_result(r, lambda, s)$r_res
})
multiple_iter_result <- sapply(r_seq, function(r) {
  single_r_result(r, lambda, s)$n_iter
})
multiple_error_result <- sapply(r_seq, function(r) {
  single_r_result(r, lambda, s)$error
})

plot(r_seq, multiple_r_result)
plot(r_seq, multiple_iter_result,
  type = "l",
  main = "lambda = 1, s = -0.905",
  ylab = "iterations to converge",
  xlab = "initial r"
)

result_0 <- pcSLOPE(1, s, r = 0.7)
pcSLOPE_target_2x2(result_0$R, result_0$D, matrix(c(1, s, s, 1), ncol = 2), lambda)
result_1 <- pcSLOPE(1, s, r = 0.5)
pcSLOPE_target_2x2(result_1$R, result_1$D, matrix(c(1, s, s, 1), ncol = 2), lambda)

pcSLOPE_target_from_r <- sapply(r_seq, function(r) {
  pcSLOPE_target_2x2(matrix(c(1, r, r, 1), ncol = 2), diag(2) / sqrt(1 + s * r), matrix(c(1, s, s, 1), ncol = 2), lambda)
})
plot(r_seq, pcSLOPE_target_from_r, type = "l")
