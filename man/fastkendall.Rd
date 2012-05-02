\name{fastkendall}
\alias{fastkendall}
\title{Kendall's tau, Goodman-Kruskal gamma}
\description{
  Fast Calculation of Kendall's tau and the Goodman-Kruskal gamma in O(n
  log(n)).  
}
\usage{
fastkendall(x, y)
}
\arguments{
  \item{x, y}{numeric vectors}
}
\details{
  Kendall's tau and the related Goodman-Kruskal gamma can be calculated
  in O(n log(n)). For large vectors, this is considerably faster than
  the O(n^2) algorithm of cor(x, y, method="kendall") in R's standard
  implementation.  Efficient algorithms are described in Knight (1966)
  and Christensen (2005), which is used here.
}
\value{
  A list of values, which may be of interest, is returned.
  \item{concordant}{}
  \item{discordant}{}
  \item{extraX}{}
  \item{extraY}{}
  \item{spare}{}
  \item{kendall.tau}{Kendall's tau}
  \item{goodman.kruskal.gamma}{Goodman-Kruskal gamma}
}
\references{
  Knight (1966)
  Christensen (2005)
}
\examples{
x = runif(1000); y=runif(1000)
fastkendall(x, y)
## for comparison
cor(x, y, method="kendall")
}