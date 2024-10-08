% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{pad_circ}
\alias{pad_circ}
\title{Pad circular data on both ends}
\usage{
pad_circ(
  data,
  circ_var,
  circ_borders = c(-90, 90),
  circ_part = 1/6,
  verbose = FALSE
)
}
\arguments{
\item{data}{data.table to pad}

\item{circ_var}{circular variable}

\item{circ_borders}{range of the circular variable}

\item{circ_part}{padding proportion}

\item{verbose}{print extra info}
}
\value{
a padded data.table
}
\description{
Pad circular data on both ends
}
\details{
Pads the data by adding a part of the data (default: 1/6th) from one end to another end. Useful to roughly account for circularity when using non-circular methods.
}
\examples{

dt <- data.table::data.table(x = runif(1000, -90, 90), y = rnorm(1000))
pad_circ(dt, "x", verbose = TRUE)

}
