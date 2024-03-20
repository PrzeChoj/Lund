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

lambdas <- 10^(seq(-2.5, 0.5, length.out = 13))



points <- matrix(ncol = 3, nrow = 0)

for(i in 1:length(lambdas)){
  print(i)
  R_out <- pcSLOPE(rep(lambdas[i], p*(p-1)/2), S)$R
  
  R_without_egdes <- R_out[-1, -1]
  should_be_zeros <- R_without_egdes[(row(R_without_egdes) > col(R_without_egdes))]
  should_be_non_zeros <- R_out[1,-1]
  for(j in 1:length(should_be_zeros)){
    points <- rbind(points, c(lambdas[i], should_be_zeros[j], 0))
  }
  for(j in 1:length(should_be_non_zeros)){
    points <- rbind(points, c(lambdas[i], should_be_non_zeros[j], 1))
  }
}

points_df <- data.frame(x = points[,1], y = points[,2], color = points[,3])

# Create scatter plot with ggplot2
ggplot(points_df, aes(x = x, y = y, color = factor(color), alpha = 0.02)) +
  geom_point() +
  scale_x_log10() +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "X (log scale)", y = "Y", title = "Scatter Plot with Color")























