remotes::install_github("PrzeChoj/LundPackage")
library(LundPackage)
#install.packages("gips")
library(ggplot2)
library(ggpubr)

source("get_matrix.R")
source("get_simulations.R")


######
# Star graph:
# In PCGLASSO article from Scandinavian Journal of Statistics, they used p = 50, n = 100, non_zero_partial_corelation = -1/sqrt(p)
p <- 50
n <- 100
non_zero_partial_corelation <- -1/sqrt(p)
star_graph <- get_matrix_star_graph(p, n, non_zero_partial_corelation)

gips:::pretty_plot_matrix(star_graph$trueK)
#gips:::pretty_plot_matrix(star_graph$trueSigma)
gips:::pretty_plot_matrix(star_graph$S)

N <- 31
lambda_min <- 0.01
lambda_max <- 1
data_star <- get_simulations(star_graph, N, lambda_min, lambda_max, SLOPE_first_multiplier = 10, progress_bar = TRUE)

jitter_width <- 0.01
dot_size <- 2
plt_LASSO <- plot_points(data_star$points_df_LASSO, non_zero_partial_corelation, "LASSO", c("black", "red"), jitter_width = jitter_width, size = dot_size)
plt_SLOPE <- plot_points(data_star$points_df_SLOPE, non_zero_partial_corelation, "SLOPE", c("blue", "purple"), jitter_width = jitter_width, size = dot_size)
plt_KL_loss <- plot_divergence(data_star$lambdas, data_star$KL_loss_LASSO, data_star$KL_loss_SLOPE, "KL loss")
plt_frob_norm <- plot_divergence(data_star$lambdas, data_star$frob_norm_LASSO, data_star$frob_norm_SLOPE, "Frob norm")

ggarrange(plt_LASSO, plt_SLOPE, plt_KL_loss, plt_frob_norm,
          ncol = 2, nrow = 2)
ggsave("./plots/star_graph.png", units = "px", width = 6000, height = 3000)


######
# cluster graph:
n <- 50
non_zero_partial_corelation <- 0.3
cluster_graph <- get_matrix_block_graph(cluster_size = 5, num_clusters = 4, corelation = non_zero_partial_corelation, n)

gips:::pretty_plot_matrix(cluster_graph$trueK)
#gips:::pretty_plot_matrix(cluster_graph$trueSigma)
gips:::pretty_plot_matrix(cluster_graph$S)

N <- 51
lambda_min <- 0.01
lambda_max <- 0.8
data_cluster <- get_simulations(cluster_graph, N, lambda_min, lambda_max, SLOPE_first_multiplier = 10, progress_bar = TRUE)

jitter_width <- 0.01
dot_size <- 2
plt_LASSO <- plot_points(data_cluster$points_df_LASSO, non_zero_partial_corelation, "PCLASSO", c("black", "red"), jitter_width = jitter_width, size = dot_size)
plt_SLOPE <- plot_points(data_cluster$points_df_SLOPE, non_zero_partial_corelation, "PCSLOPE", c("blue", "purple"), jitter_width = jitter_width, size = dot_size)
plt_KL_loss <- plot_divergence(data_cluster$lambdas, data_cluster$KL_loss_LASSO, data_cluster$KL_loss_SLOPE, "KL loss")
plt_frob_norm <- plot_divergence(data_cluster$lambdas, data_cluster$frob_norm_LASSO, data_cluster$frob_norm_SLOPE, "Frob norm")

ggarrange(plt_LASSO, plt_SLOPE, plt_KL_loss, plt_frob_norm,
          ncol = 2, nrow = 2)
ggsave("./plots/cluster_graph.png", units = "px", width = 6000, height = 3000)


