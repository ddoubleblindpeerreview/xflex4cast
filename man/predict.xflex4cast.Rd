% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flexForecast_model.R, R/xflex4cast_model.R
\name{predict.xflex4cast}
\alias{predict.xflex4cast}
\title{Predict from xflex4cast Model}
\usage{
\method{predict}{xflex4cast}(model, newdata, indexes = NULL)

\method{predict}{xflex4cast}(model, newdata, indexes = NULL)
}
\arguments{
\item{model}{A fitted \code{xflex4cast} model.}

\item{newdata}{A \code{data.frame} or \code{data.table} of new data to predict on.}

\item{indexes}{Names of identifying columns to track predictions. Default: \code{index} with row numbers.}
}
\value{
A \code{data.table} in long format with columns: indexes, quantile, and value.

A \code{data.table} in long format with columns: indexes, quantile, and value.
}
\description{
Predict quantiles from a fitted \code{xflex4cast} model. Returns long-format table.

Predict quantiles from a fitted \code{xflex4cast} model. Returns long-format table.
}
