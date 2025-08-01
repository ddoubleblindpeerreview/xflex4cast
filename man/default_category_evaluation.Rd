% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation_metrics.R
\name{default_category_evaluation}
\alias{default_category_evaluation}
\title{Categorical forecast evaluation for P4R models}
\usage{
default_category_evaluation(
  predictions,
  data,
  categories,
  indexes = list(valid_time = "dtm", base_time = "dtm_basetime"),
  lead_time_units = "hours",
  ...
)
}
\arguments{
\item{predictions}{probabilistic fault forecasts, quantiles in long format
\code{data.table}}

\item{data}{actual fault data with date-time column(s) and column
\code{actual} containing outcome corresponding to predictions. A merge will be
attempted between \code{predictions} and \code{data} on all column names in both
\code{indexes} and \code{data}.}

\item{categories}{category definitions, a vector of length 2
\code{c(green-amber,amber-red)}.}

\item{indexes}{a list with elements \code{valid_time} and \code{base_time}.
\code{base_time} my be \code{NULL}.}

\item{lead_time_units}{passed to \code{difftime()} to specify units for forecast
lead-time in results}

\item{...}{Additional arguments passed to \code{get_confusion_matrix()}}
}
\value{
A list containing confusion matrices and \code{pROC} objects.
}
\description{
Categorical forecast evaluation for P4R models
}
\details{
Calculates confusion matrix by regulatory year and lead-time, and
ROC by lead-time.
}
