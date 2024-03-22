#' star_graph <- get_matrix_star_graph(4, 60)
#' star_graph$S
#' star_graph$trueSigma
#' star_graph$trueK
get_matrix_star_graph <- function(p, n, non_zero_partial_corelation){
  trueK <- diag(p)
  trueK[1, -1] <- non_zero_partial_corelation
  trueK[-1, 1] <- non_zero_partial_corelation
  
  trueSigma <- solve(trueK)
  
  Z <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = trueSigma)
  S <- (t(Z) %*% Z) / n
  
  list("trueSigma" = trueSigma, "trueK" = trueK, "S" = S)
}

#' star_graph <- get_matrix_block_graph(4, 3, 30)
#' star_graph$S
#' star_graph$trueSigma
#' star_graph$trueK
get_matrix_block_graph <- function(cluster_size, num_clusters, corelation, n){
  p <- cluster_size * num_clusters
    
  one_cluster <- diag(cluster_size)
  one_cluster[col(one_cluster) != row(one_cluster)] <- corelation
  big_matrix <- diag(num_clusters)
  
  trueK <- kronecker(big_matrix, one_cluster)
  
  trueSigma <- solve(trueK)
  
  Z <- MASS::mvrnorm(n = n, mu = numeric(p), Sigma = trueSigma)
  S <- (t(Z) %*% Z) / n
  
  list("trueSigma" = trueSigma, "trueK" = trueK, "S" = S)
}
