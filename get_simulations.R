single_update_points <- function(S, points, lambdas, lambdas_i, places_should_be_zero, places_should_be_non_zero) {
  my_out <- pcSLOPE(lambdas, S)
  R_out <- my_out$R
  
  should_be_zeros <- R_out[places_should_be_zero]
  should_be_non_zeros <- R_out[places_should_be_non_zero]
  
  for(j in 1:length(should_be_zeros)){
    points <- rbind(points, c(lambdas_i, should_be_zeros[j], 0))
  }
  for(j in 1:length(should_be_non_zeros)){
    points <- rbind(points, c(lambdas_i, should_be_non_zeros[j], 1))
  }
  
  attr(points, "R") <- R_out
  
  points
}

from_norm <- function(A) {
  sqrt(sum(A^2))
}

get_simulations <- function(graph, N, lambda_min, lambda_max,
                            SLOPE_first_multiplier = 1, progress_bar = TRUE) {
  S <- graph$S
  trueSigma <- graph$trueSigma
  trueK <- graph$trueK
  
  p <- dim(S)[1]
  
  possible_values <- as.numeric(names(table(trueK[row(trueK) != col(trueK)])))
  possible_values <- setdiff(possible_values, 0) # we know there should be zero
  
  if(length(possible_values) > 1){
    rlang::abort("The code is assuming there is only one non-zero off-diagonal in true K")
  }
  
  places_should_be_zero <- (abs(trueK) < 0.000001)
  places_should_be_non_zero <- (abs(trueK - possible_values) < 0.000001)
  
  lambdas <- 10^(seq(log(lambda_min, 10), log(lambda_max, 10), length.out = N))
  
  KL_loss_LASSO <- numeric(N)
  KL_loss_SLOPE <- numeric(N)
  
  frob_norm_LASSO <- numeric(N)
  frob_norm_SLOPE <- numeric(N)
  
  points_LASSO <- matrix(ncol = 3, nrow = 0)
  points_SLOPE <- matrix(ncol = 3, nrow = 0)
  
  if (progress_bar) {
    pb <- txtProgressBar(min = 0, max = length(lambdas) * 2, initial = 0) 
  }
  
  for(i in 1:length(lambdas)){
    if (progress_bar) {
      setTxtProgressBar(pb, i * 2)
    }
    
    # LASSO:
    lambdas_LASSO <- rep(lambdas[i], p*(p-1)/2)
    points_LASSO <- single_update_points(S, points_LASSO, lambdas_LASSO, lambdas[i], places_should_be_zero, places_should_be_non_zero)
    
    KL_loss_LASSO[i] <- KL_loss(trueSigma, trueK, attr(points_LASSO, "R"))
    frob_norm_LASSO[i] <- from_norm(attr(points_LASSO, "R") - trueK)
    
    if (progress_bar) {
      setTxtProgressBar(pb, i * 2 + 2)
    }
    
    # SLOPE:
    lambdas_SLOPE <- lambdaBH(p, n, alpha = 0.000000001, l_max = lambdas[i])
    #lambdas_SLOPE <- lambdaLinear(p, l_max = 2*lambdas[i], l_min = 0.2*lambdas[i])
    lambdas_SLOPE[1] <- lambdas_SLOPE[1] * SLOPE_first_multiplier
    #lambdas_SLOPE[2] <- lambdas_SLOPE[2] * 3
    #lambdas_SLOPE[3] <- lambdas_SLOPE[3] * 1.5
    
    points_SLOPE <- single_update_points(S, points_SLOPE, lambdas_SLOPE, lambdas[i], places_should_be_zero, places_should_be_non_zero)
    
    KL_loss_SLOPE[i] <- KL_loss(trueSigma, trueK, attr(points_SLOPE, "R"))
    frob_norm_SLOPE[i] <- from_norm(attr(points_SLOPE, "R") - trueK)
  }
  
  points_df_LASSO <- data.frame(x = points_LASSO[,1], y = points_LASSO[,2], color = points_LASSO[,3])
  points_df_SLOPE <- data.frame(x = points_SLOPE[,1], y = points_SLOPE[,2], color = points_SLOPE[,3])
  
  list("points_df_LASSO" = points_df_LASSO, "points_df_SLOPE" = points_df_SLOPE,
       "KL_loss_LASSO" = KL_loss_LASSO, "KL_loss_SLOPE" = KL_loss_SLOPE,
       "frob_norm_LASSO" = frob_norm_LASSO, "frob_norm_SLOPE" = frob_norm_SLOPE,
       "lambdas" = lambdas)
}

plot_points <- function(data_points, non_zero_value, title, my_colors = c("black", "red"), jitter_width = 0, ylim_top = NULL, size = 2){
  if (is.null(ylim_top)) {
    ylim_top <- max(abs(data_points$y)) * 1.1
  }
  
  ggplot(data_points, aes(x = x, y = y, color = factor(color))) +
    geom_jitter(width = jitter_width, size = size, alpha = 0.3) +
    scale_x_log10(n.breaks = 6) +
    ylim(-ylim_top, ylim_top) +
    scale_color_manual(values = my_colors, name = "Should be", labels = c("0", sprintf("%.2f", non_zero_value))) +
    labs(x = "Regularization parameter", y = "Partial correlation estimates", title = title) +
    geom_hline(yintercept = 0, linetype="dashed", color = my_colors[1]) +
    geom_hline(yintercept = non_zero_value, linetype="dashed", color = my_colors[2])
}

plot_divergence <- function(lambdas, value_LASSO, value_SLOPE, y_axis_name){
  data <- data.frame(x = rep(lambdas, 2),
                     y = c(value_LASSO, value_SLOPE),
                     group = rep(c("LASSO", "SLOPE"), each = length(lambdas)))
  
  max_y <- max(max(value_LASSO), max(value_SLOPE))
  
  # Plot using ggplot2
  ggplot(data, aes(x = x, y = y, color = group)) +
    geom_point() +
    geom_line() +
    scale_x_log10(n.breaks = 6) +
    ylim(0, max_y) +
    labs(x = "Regularization parameter", y = y_axis_name) +
    scale_color_manual(values = c("red", "purple"))
}
