#include "rcpp_comembership.h"

SEXP rcpp_comembership(SEXP labels) {
  using namespace Rcpp;
    
  Rcpp::NumericVector cluster_labels(labels);
  int n = cluster_labels.size();

  // The comembership vector is of length "n choose 2".
  Rcpp::NumericVector comembership(n * (n-1) / 2);

  // The comembership index.
  int idx_comembership = 0;

  // Traverse through all pairs of observations to identify comemberships.
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      if (cluster_labels[i] == cluster_labels[j]) {
        comembership[idx_comembership] = 1;
      }
      idx_comembership++;
    }
  }
  
  return comembership;
}
