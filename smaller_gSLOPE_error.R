gslope_standardized_new <- function(
    sample_cov,
    lambdas,
    rho = 1.1,
    max_iter = 500,
    epsilon = 1e-04,
    verbose = FALSE) {
  if (!(nrow(sample_cov) == ncol(sample_cov))) {
    stop("Covariance matrix must be square.")
  }

  # Parameters initialization
  Z <- sample_cov # Initialized to zero, probably it is the best choice
  Y <- Z
  diag(Y) <- 1
  X <- diag(nrow(sample_cov))

  for (iter in 1:max_iter) {
    if (verbose) {
      print(diag(X))
    }

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
    if (primal_residual < epsilon && dual_residual < epsilon) {
      break
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






inputMatrix <- structure(
  c(
    2.72851295826233, -2.46930422608904,
    -2.46930422608904, 2.7285129557466
  ),
  dim = c(2L, 2L)
)
lambda <- 1
gslope_standardized_new(inputMatrix, lambda, verbose = TRUE)
