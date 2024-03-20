remotes::install_github("PrzeChoj/LundPackage")
library(LundPackage)
library(ggplot2)

# star graph example
p <- 50
n <- 100
trueK <- diag(p)
trueK[1, -1] <- -1/sqrt(p)
trueK[-1, 1] <- -1/sqrt(p)

trueSigma <- solve(trueK)

Z <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = trueSigma)
S <- (t(Z) %*% Z) / n

#install.packages("gips")
gips:::pretty_plot_matrix(trueSigma)
gips:::pretty_plot_matrix(S)

N <- 5
lambdas <- 10^(seq(-0.5, 0.5, length.out = N))
KL_loss <- numeric(N)
n_iters <- numeric(N)


points <- matrix(ncol = 3, nrow = 0)

for(i in 1:length(lambdas)){ #NOTE: For max_iter = 100, it does not converge sometimes yet
  print(i)
  lambdas_LASSO <- rep(lambdas[i], p*(p-1)/2)
  lambdas_SLOPE <- lambdaBH(p, n, alpha = 0.000000001, l_max = lambdas[i])
  #lambdas_SLOPE <- lambdaLinear(p, l_max = lambdas[i])
  my_out <- pcSLOPE(lambdas_LASSO, S)
  R_out <- my_out$R
  n_iters[i] <- my_out$iter
  
  R_without_egdes <- R_out[-1, -1]
  should_be_zeros <- R_without_egdes[(row(R_without_egdes) > col(R_without_egdes))]
  should_be_non_zeros <- R_out[1,-1]
  for(j in 1:length(should_be_zeros)){
    points <- rbind(points, c(lambdas[i], should_be_zeros[j], 0))
  }
  for(j in 1:length(should_be_non_zeros)){
    points <- rbind(points, c(lambdas[i], should_be_non_zeros[j], 1))
  }
  
  # SLOPE:
  lambdas_SLOPE[1] <- lambdas_SLOPE[1] * 10
  my_out <- pcSLOPE(lambdas_SLOPE, S)
  R_out <- my_out$R
  n_iters[i] <- my_out$iter
  
  R_without_egdes <- R_out[-1, -1]
  should_be_zeros <- R_without_egdes[(row(R_without_egdes) > col(R_without_egdes))]
  should_be_non_zeros <- R_out[1,-1]
  for(j in 1:length(should_be_zeros)){
    points <- rbind(points, c(lambdas[i], should_be_zeros[j], 2))
  }
  for(j in 1:length(should_be_non_zeros)){
    points <- rbind(points, c(lambdas[i], should_be_non_zeros[j], 3))
  }
  
  KL_loss[i] <- KL_loss(trueSigma, trueK, my_out$K)
}

points_df <- data.frame(x = points[,1], y = points[,2], color = points[,3])
ggplot(points_df, aes(x = x, y = y, color = factor(color), alpha = 0.02, size =2)) +
  geom_point() +
  scale_x_log10() +
  ylim(-0.5, 0.5) +
  scale_color_manual(values = c("black", "red", "blue", "purple")) +
  labs(x = "Regularization parameter", y = "Partial correlation estimates", title = "SLOPE")


data <- data.frame(ox = lambdas, oy = KL_loss)
ggplot(data, aes(x = ox, y = oy)) +
  geom_point() +
  #scale_x_log10() +
  labs(x = "Regularization parameter", y = "KL loss",
       title = "PCGLASSO")

data <- data.frame(ox = lambdas, oy = n_iters)
ggplot(data, aes(x = ox, y = oy)) +
  geom_point() +
  #scale_x_log10() +
  labs(x = "Regularization parameter", y = "number of iteration for SLOPE solver",
       title = "PCGLASSO")





