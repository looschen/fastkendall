\name{fastkendall}
\alias{fastkendall}
\title{Kendall's \eqn{\tau}, Goodman-Kruskal \eqn{\gamma}}
\description{
  Fast Calculation of Kendall's \eqn{\tau} and the Goodman-Kruskal
  \eqn{\gamma} in \eqn{O(n log(n))}.
}
\usage{
fastkendall(x, y)
}
\arguments{
  \item{x, y}{numeric vectors}
}
\details{
  Kendall's \eqn{\tau} and the related Goodman-Kruskal \eqn{\gamma} are
  rank correlaton coefficients, which compare pairs \code{(x[i], x[j])}
  and \code{(y[i], y[j])} for all \code{i} and \code{j}.

  These rank correlation coeffients can be calculated in \eqn{O(n
  log(n))}. For large vectors, this is considerably faster than the
  \eqn{O(n^2)} algorithm of \code{cor(x, y, method="kendall")} in R's
  standard implementation.  Efficient algorithms are described in
  \cite{Knight (1966)} and \cite{Christensen (2005)}. The approach of
  the latter is used here.
}
\value{
  A list of values of interest is returned.
  \item{concordant}{number of pairs with \code{x[i] < x[j] && y[i] < y[j]}\cr (or inverted "<")}
  \item{discordant}{number of pairs with \code{x[i] < x[j] && y[i] > y[j]}\cr (or inverted "<", ">")}
  \item{extra.x}{number of pairs with \code{x[i] != x[j] && y[i] == y[j]}}
  \item{extra.y}{number of pairs with \code{x[i] == x[j] && y[i] != y[j]}}
  \item{spare}{number of pairs with \code{x[i] == x[j] && y[i] == y[j]}}
  \item{kendall.tau}{Kendall's \eqn{\tau},\cr
    \preformatted{tau = (concordant - discordant) / (
      sqrt(concordant + discordant + extra.x) *
      sqrt(concordant + discordant + extra.y))}}
  \item{goodman.kruskal.gamma}{Goodman-Kruskal \eqn{\gamma},\cr
    \preformatted{gamma = (concordant - discordant) /
      (concordant + discordant)}
  }
}
\author{
  Karsten Looschen, \email{karsten.looschen@gmail.com}
}
\references{
  % Knight (1966)
  W. R. Knight. A Computer Method for Calculating Kendall's Tau with
  UngroupedData. Journal of the American Statistical Association,
  61(314):436--439, 1966.
  
  % Christensen (2005)
  D. Christensen. Fast algorithms for the calculation of Kendall's \eqn{\tau}
  . Computational Statistics, 20:51--62, 2005.
}
\seealso{
  \code{\link{cor}}
}
\examples{
x = runif(1000); y=runif(1000)
fastkendall(x, y)
## compare with R's implementation
cor(x, y, method="kendall")

## correlate data with ties (probably)
x = sample(10, 100, replace=TRUE); y = sample(10, 100, replace=TRUE)
fastkendall(x, y)
}
\keyword{univar}