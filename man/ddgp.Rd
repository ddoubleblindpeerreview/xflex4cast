% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tail_computations.R
\name{ddgp}
\alias{ddgp}
\title{Probability Mass Function of the Discrete Generalized Pareto Distribution}
\usage{
ddgp(x, shape, scale)
}
\arguments{
\item{x}{Numeric vector of non-negative integers at which to evaluate the PMF.}

\item{shape}{Numeric scalar. Shape parameter of the DGP distribution.}

\item{scale}{Numeric scalar. Scale parameter of the DGP distribution.}
}
\value{
Numeric vector of PMF values corresponding to \code{x}.
}
\description{
Computes the probability mass function (PMF) values of the Discrete Generalized Pareto (DGP) distribution at given points.
}
\details{
The DGP distribution is a discrete analogue of the continuous Generalized Pareto distribution. When the shape parameter is near zero, the PMF simplifies to differences of exponential terms; otherwise, it follows the standard DGP formula.

The support of the distribution is \eqn{x \in \{0,1,2,...\}}.
}
\examples{
ddgp(x = 0:5, shape = 0.3, scale = 2)

}
