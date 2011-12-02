library('testthat')
library('clusteval')

set.seed(42)

iris_x <- data.matrix(iris[, -5])
iris_y <- iris[, 5]

# TODO: Create consensus plots for the various clustering algorithms and various K
# TODO: Create a joint consensus clustering algorithm for various K
# TODO: Make unit test from this example
consensus_out <- consensus(x = iris_x, num_clusters = 3, method = "kmeans",
                           num_resamples = 100, subsample_size = 50)
print(plot_consensus(consensus_out$consensus_mat))
