% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{circ_mean_rad}
\alias{circ_mean_rad}
\alias{circ_mean_180}
\alias{circ_mean_360}
\title{Circular mean}
\usage{
circ_mean_rad(x, na.rm = FALSE)

circ_mean_180(x, na.rm = FALSE)

circ_mean_360(x, na.rm = FALSE)
}
\arguments{
\item{x}{vector of values}

\item{na.rm}{a logical value indicating whether NA values should be removed before the computation proceeds}
}
\value{
mean of values in the vector
}
\description{
Circular mean
}
\section{Functions}{
\itemize{
\item \code{circ_mean_rad()}: circular mean in 2pi space

\item \code{circ_mean_180()}: circular mean in 180° space (e.g., line orientation)

\item \code{circ_mean_360()}: circular mean in 360° space

}}
\examples{
x <- runif(1000, -pi, pi)
mean(x)
circ_mean_rad(x)

}
