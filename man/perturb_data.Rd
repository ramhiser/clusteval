\name{perturb_data}
\alias{perturb_data}
\title{Perturbs a multivariate data set stored in a matrix.}
\usage{
  perturb_data(x, noise = 0.1)
}
\arguments{
  \item{x}{TODO}

  \item{noise}{TODO}
}
\value{
  TODO
}
\description{
  On the blog "X-Informatics" there is a post called "Add
  noise to data" where a method is described to add noise
  to a data set. The method finds the smallest dimension of
  the bounding box enclosing the data. Take the square root
  of this, multiply by 'noise' (which is defaulted to
  0.02), and use this as the standard deviation when
  generating a normal random variate.
}
\details{
  The idea originates in the paper: Balasubramanian, M.,
  Schwartz, E.L., Tenenbaum, J.B., de Silva, V. & Langford,
  J.C. "The Isomap Algorithm and Topological Stability".
  Science 295, 7a (2002).

  The blog post is located:
  http://c13s.wordpress.com/2010/04/13/add-noise-to-data/

  A note about the minimum bounding box idea is located:
  http://en.wikipedia.org/wiki/Minimum_bounding_box
}

