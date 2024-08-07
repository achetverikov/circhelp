% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{predict.circ_loess}
\alias{predict.circ_loess}
\title{Compute predictions for circular LOESS}
\usage{
\method{predict}{circ_loess}(object, newdata, ...)
}
\arguments{
\item{object}{a circular LOESS object}

\item{newdata}{a data.frame with a variable x on which the predictions are computed}

\item{...}{other arguments (passed to circ_loess)}
}
\value{
a data.frame with predictions
}
\description{
Compute predictions for circular LOESS
}
\keyword{internal}
