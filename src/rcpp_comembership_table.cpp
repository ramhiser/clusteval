#include "rcpp_comembership_table.h"

SEXP rcpp_comembership_table(SEXP labels1, SEXP labels2) {
  using namespace Rcpp;
    
  Rcpp::NumericVector cluster_labels1(labels1);
  Rcpp::NumericVector cluster_labels2(labels2);

  int n1 = cluster_labels1.size(), n2 = cluster_labels2.size();
  if (n1 != n2) {
    throw std::range_error("The two vectors of cluster labels must be of equal length.");
  }

  // The counts of comembership pairs.
  // n_11: the number of comemberships in both partitions
  // n_10: the number of comemberships in clustering 1 but not in clustering 2
  // n_01: the number of comemberships in clustering 2 but not in clustering 1
  // n_00: the number of non-comemberships in both partitions
  int n_11 = 0, n_10 = 0, n_01 = 0, n_00 = 0;

  // Flags that indicate if the current pair of cluster labels in clusterings 1
  // and 2 if are comemberships.
  bool comembership1 = false, comembership2 = false;

  // Traverse through all pairs of observations to identify comemberships.
  for (int i = 0; i < n1; i++) {
    for (int j = i + 1; j < n1; j++) {
      // If either of the clusterings have a comembership, set the comembership
      // flag as true.
      comembership1 = (cluster_labels1[i] == cluster_labels1[j]);
      comembership2 = (cluster_labels2[i] == cluster_labels2[j]);

      if (comembership1 && comembership2) {
        n_11++;
      } else if (comembership1 && !comembership2) {
        n_10++;
      } else if (!comembership1 && comembership2) {
        n_01++;
      } else { // if (!comembership1 && !comembership2)
        n_00++;
      }
    }
  }
  
  // Returns a list that contains the 2x2 contingecy table results.
  return Rcpp::List::create(Rcpp::Named("n_11") = n_11,
                            Rcpp::Named("n_10") = n_10,
                            Rcpp::Named("n_01") = n_01,
                            Rcpp::Named("n_00") = n_00);
}
