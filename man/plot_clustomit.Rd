\name{plot_clustomit}
\alias{plot_clustomit}
\title{Plot of clustomit Results}
\usage{
  plot_clustomit(obj, binwidth = NULL)
}
\arguments{
  \item{obj}{object of class 'clustomit' with results from
  clustomit}

  \item{binwidth}{the histogram bin width for the resulting
  ggplot object. Use ggplot2's defaults if NULL.}
}
\value{
  a ggplot2 graphical summary of the results
}
\description{
  TODO
}
\details{
  Example: B <- 500 num_clusters <- c(3, 4) cluster_methods
  <- c("kmeans") similarity_methods <- c("rand", "jaccard")
  grid <- expand.grid(K = num_clusters, cluster =
  cluster_methods, similarity = similarity_methods)

  out <- mapply(clustomit_landon, K = grid$K,
  cluster_method = grid$cluster, similarity_method =
  grid$similarity, MoreArgs = list(x = iris_x, B = B),
  SIMPLIFY = FALSE)

  binwidth <- NULL obj <- out[[4]] plot_clustomit(obj)
}

