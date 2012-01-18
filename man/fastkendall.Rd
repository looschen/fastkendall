\name{fastkendall}
\alias{fastkendall}
\title{Fast Calculation of Kendall's tau and the Goodman-Kruskal-gamma}
\description{
  Calculate Stuff much faster than cor(...,
  method="kendall")
}
\usage{
fastkendall(x, y)
}
\arguments{
  \item{x, y} numeric vectors
}

\examples{
x = rnorm(1e4); y=rnorm(1e4)
fk = fastkendall(x, y)
stopifnot(fk$kendall.tau == cor(x, y, method="kendall"))
}