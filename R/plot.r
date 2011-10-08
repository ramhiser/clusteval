#' Plot of clustomit Results
#'
#' TODO
#'
#' Example:
#' B <- 500
#' num_clusters <- c(3, 4)
#' cluster_methods <- c("kmeans")
#' similarity_methods <- c("rand", "jaccard")
#' grid <- expand.grid(K = num_clusters, cluster = cluster_methods, similarity = similarity_methods)
#' 
#' out <- mapply(clustomit_landon, K = grid$K, cluster_method = grid$cluster, similarity_method = grid$similarity, MoreArgs = list(x = iris_x, B = B), SIMPLIFY = FALSE)
#' 
#' binwidth <- NULL
#' obj <- out[[4]]
#' plot_clustomit(obj)
#'
#' @param obj object of class 'clustomit' with results from clustomit
#' @param binwidth the histogram bin width for the resulting ggplot object. Use ggplot2's defaults if NULL.
#' @return a ggplot2 graphical summary of the results
plot_clustomit <- function(obj, binwidth = NULL) {
  if(!inherits(obj, "clustomit")) {
    stop("obj not of class 'clustomit'")
  }
  simil_string <- switch(obj$similarity_method,
    "rand" = "Rand",
    "jaccard" = "Jaccard",
    "adjusted_rand" = "Adjusted Rand"
  )

  clust_string <- switch(obj$cluster_method,
    "kmeans" = "K-Means",
    "singlelinkage" = "Hierarchical (Single-Linkage)",
    "pam" = "PAM",
    "diana" = "Diana",
    "mclust" = "Model-based"
  )

  plot_title <- "Bootstrapped Sampling Distribution\n"
  plot_title <- paste(plot_title, clust_string, "Clustering with K =", obj$K)
  
  plot_strings <- list()
  plot_strings$obsvalue <- paste("Observed Value:", round(obj$observed_value, 3))
  obsvalue_string <- paste("Observed Value:", round(obj$observed_value, 3))
  
  scores <- data.frame(score = obj$scores)

  p <- ggplot(scores, aes(x = score))
  p <- p + geom_histogram(color = "darkgreen", fill = "white", binwidth = binwidth)
  p <- p + geom_vline(xintercept = obj$observed_value, color = "darkgreen", size = 1)
  #p <- p + geom_text(aes(x = obj$observed_value, y = 0, label = plot_strings$obsvalue), color = "darkgreen", size = 4, hjust = -0.1, vjust = 2)
  p <- p + opts(title = plot_title) + xlab(paste(simil_string, "Similarity"))
  p
}