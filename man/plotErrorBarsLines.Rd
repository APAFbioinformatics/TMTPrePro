\name{plotErrorBarsLines}
\alias{plotErrorBarsLines}
\title{
Utility for clustering plots to plot lines and an overall trend
}
\description{
Prints faint lines for each profile, and a mean/error bars
}
\usage{
plotErrorBarsLines(v, barSizes, lines, labels = NULL, col = "blue", ylim, ...)
}
\arguments{
  \item{v}{
	Overall trend, to be printed solid, length n
}
  \item{barSizes}{
	Size of the error bars, length n
}
  \item{lines}{
	Matrix of n columns, and as many rows as lines
}
  \item{labels}{
	Labels to be printed on the x axis, length n
}
  \item{col}{
	Colour for main trend line
}
  \item{ylim}{
	Can specfy limits so several graphs are on the same scale
}
  \item{\dots}{
	Additional parameters to pass in
}
}
\details{
}
\value{
No returned value; plot only.
}
\references{
}
\author{
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{help}}, ~~~
}
\examples{
mat = matrix(rnorm(100), 10)
plotErrorBarsLines(apply(mat,1,FUN=mean), apply(mat,1,FUN=sd), 
	lines=mat, col="red", main="A random plot", xlab="Some label")


}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
