# Wrapper for the Folkes-Mallows index for comparing two clusterings of the same data set.
folkes_mallows_wrapper <- function(cl1, cl2) {
	clv.Folkes.Mallows(std.ext(cl1, cl2))
}

# Wrapper for the Phi index for comparing two clusterings of the same data set.
phi_wrapper <- function(cl1, cl2) {
	clv.Phi(std.ext(cl1, cl2))
}

# Wrapper for the Russel-Rao index for comparing two clusterings of the same data set.
russel_rao_wrapper <- function(cl1, cl2) {
	clv.Russel.Rao(std.ext(cl1, cl2))
}