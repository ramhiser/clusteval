\name{intraclass_cov}
\alias{intraclass_cov}
\title{Construct an intraclass covariance matrix.}
\usage{
  intraclass_cov(p, rho, sigma2 = 1)
}
\arguments{
  \item{p}{the dimension of the matrix}

  \item{rho}{the intraclass covariance (correlation)
  constant}

  \item{sigma2}{the coefficient of the intraclass
  covariance matrix}
}
\value{
  an intraclass covariance matrix matrix of size p
}
\description{
  We define a \eqn{p \times p} intraclass covariance
  (correlation) matrix to be \eqn{\Sigma_m = \sigma^2 (1 -
  \rho) J_p + \rho I_p}, where \eqn{-(p-1)^{-1} < \rho <
  1}, \eqn{I_p} is the \eqn{p \times p} identity matrix,
  and \eqn{J_p} denotes the \eqn{p \times p} matrix of
  ones.
}

