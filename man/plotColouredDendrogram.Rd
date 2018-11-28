\name{plotColouredDendrogram}
\alias{plotColouredDendrogram}
\title{
HClust utility to plot coloured end leaves on a dendrogram
}
\description{
}
\usage{
plotColouredDendrogram(result, glabel, clabel = NULL, ...)
}
\arguments{
  \item{result}{
	Result returned by hclust
}
  \item{glabel}{
    Sample label
}
  \item{clabel}{
    Colour groups, typically resulting from cutting the clustering tree somewhere
}
  \item{\dots}{
  Other parameters passed to plot, such as plot title etc
}
}
\details{
}
\value{
Plotting only, no returned result
}
\references{
}
\author{
}
\note{
}

\seealso{
\code{\link{hclust}}
}
\examples{
result = hclust(dist(log(iris[,-5])))
plotColouredDendrogram(result, glabel=iris[,5], clabel=iris[,5])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
