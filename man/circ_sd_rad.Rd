% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{circ_sd_rad}
\alias{circ_sd_rad}
\alias{circ_sd_360}
\alias{circ_sd_180}
\title{Circular standard deviation}
\usage{
circ_sd_rad(x, na.rm = FALSE)

circ_sd_360(x, na.rm = FALSE)

circ_sd_180(x, na.rm = FALSE)
}
\arguments{
\item{x}{vector of angles}

\item{na.rm}{a logical value indicating whether NA values should be removed before the computation proceeds}
}
\value{
standard deviation of values in the vector
}
\description{
Circular standard deviation
}
\section{Functions}{
\itemize{
\item \code{circ_sd_rad()}: SD of angles in radians

\item \code{circ_sd_360()}: SD of angles in 360 degree space

\item \code{circ_sd_180()}: SD of angles in 180 degree space

}}
\examples{
circ_sd_rad(rnorm(50))
circ_sd_180(rnorm(50))

}
