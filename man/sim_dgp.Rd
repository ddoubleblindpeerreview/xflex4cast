% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_functions.R
\name{sim_dgp}
\alias{sim_dgp}
\title{Simulate Data from a Discrete Generalized Pareto Process}
\usage{
sim_dgp(n = 100, xi = 0.3)
}
\arguments{
\item{n}{Integer. Number of samples to generate.}

\item{xi}{Numeric. Shape parameter for the DGP distribution.}
}
\value{
A list with:
\describe{
\item{data}{A \code{data.table} with covariates \code{x1}, \code{x2}, and simulated outcome \code{y}.}
\item{scale_parameter_dgp}{Vector of scale parameters for the DGP model per observation.}
}
}
\description{
Generates data based solely on a DGP model where the scale parameter depends on covariates. Useful as a pure heavy-tail generator.
}
\details{
The covariates \code{x1} and \code{x2} are generated from scaled beta distributions and used to model the scale of the DGP distribution via a log-linear formulation.
}
\examples{
out <- sim_dgp()
plot(table(out$data$y), xlab = "y", ylab = "Frequency")

}
