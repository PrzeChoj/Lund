remotes::install_github("PrzeChoj/LundPackage")
library(LundPackage)

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


N_l <- 10
l_seq <- seq(0, 3 / 2, length.out = N_l + 2)[2:(N_l + 1)]
N_s <- 20
s_seq <- seq(-1, 1, length.out = N_s + 2)[2:(N_s + 1)]

args_grid <- tidyr::crossing(lambda = l_seq, s = s_seq)

single_r_result <- function(r, lambda, s) {
  pcSLOPE_res <- pcSLOPE_2x2(lambda, s, r)
  list(r_res = pcSLOPE_res$R[1, 2], n_iter = pcSLOPE_res$iter, error = pcSLOPE_res$error)
}

r_gap <- function(lambda, s) {
  N_r <- 10
  r_seq <- seq(-1, 1, length.out = N_r + 2)[2:(N_r + 1)]
  
  multiple_r_result <- sapply(r_seq, function(r) {
    single_r_result(r, lambda, s)$r_res
  })
  
  (max(multiple_r_result) - min(multiple_r_result))
}

start_time <- Sys.time()
result <- purrr::pmap_dbl(args_grid, r_gap) # 15 seconds
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

result_0 <- LundPackage::pcSLOPE_2x2(1, s, r = 0.7)
pcSLOPE_target_2x2(result_0$R, result_0$D, matrix(c(1, s, s, 1), ncol = 2), lambda)
result_1 <- pcSLOPE_2x2(1, s, r = 0.5)
pcSLOPE_target_2x2(result_1$R, result_1$D, matrix(c(1, s, s, 1), ncol = 2), lambda)

pcSLOPE_target_from_r <- sapply(r_seq, function(r) {
  pcSLOPE_target_2x2(matrix(c(1, r, r, 1), ncol = 2), diag(2) / sqrt(1 + s * r), matrix(c(1, s, s, 1), ncol = 2), lambda)
})
plot(r_seq, pcSLOPE_target_from_r, type = "l")
