library('bootclust')

context("Three bivariate Gaussian with equal, identity cov matrices")

normal <- function(N = 30, p = 2, delta = 1, probs = c(1,1,1)/3) {
 stopifnot(sum(probs) == 1) 
}

k <- 3
delta <- 1
p <- 2 
n1 <- 20
n2 <- 20
n3 <- 20
mu1 <- 1 * rep(1,p) * delta
mu2 <- 2 * rep(1,p) * delta
mu3 <- 3 * rep(1,p) * delta
Sigma <- diag(p)

x <- rbind(rmvnorm(n1, mu1, Sigma), rmvnorm(n2, mu2, Sigma), rmvnorm(n3, mu3, Sigma))

# Distances for clustering
euc_dist_x <- dist(x)
man_dist_x <- dist(x, method = "manhattan")

hc_ave_euc <- cutree(hclust(euc_dist_x, "ave"), k = 3)
hc_ave_man <- cutree(hclust(man_dist_x, "ave"), k = 3)
hc_single_euc <- cutree(hclust(euc_dist_x, "single"), k = 3)
hc_single_man <- cutree(hclust(man_dist_x, "single"), k = 3)
hc_compl_euc <- cutree(hclust(euc_dist_x, "compl"), k = 3)
hc_compl_man <- cutree(hclust(man_dist_x, "compl"), k = 3)
kmeans_out <- kmeans(x, centers = 3)
