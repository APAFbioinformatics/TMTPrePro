\name{HClust}
\alias{HClust}
\title{
Function to do hierarchical clustering
}
\description{
Hierarchical clustering using hclust, with added coloured end dendrogram.
}
\usage{
HClust(data = data, metric = c("euclidean", "manhattan", "pearsonCorrelation"), scale = FALSE, method = c("single", "complete", "average"), basefile = "hclust", glabel = row.names(data), clabel = NULL, cutNumber = NULL, cutHeight = NULL, plot = TRUE)
}
\arguments{
  \item{data}{
	Matrix to cluster
}
  \item{metric}{
	A choice of euclidean, manhattan or correlation
}
  \item{scale}{
	Should the matrix be scaled by row first
}
  \item{method}{
	A choice of single, complete or average
}
  \item{basefile}{
	Name of resulting plot
}
  \item{glabel}{
	Labels for the data rows
}
  \item{clabel}{
	Labels for the data column
}
  \item{cutNumber}{
	In how many pieces to cut the resulting dendrogram
}
  \item{cutHeight}{
	At what height to cut
}
  \item{plot}{
	TRUE/FALSE, should dendrogram be generated
}
}
\details{
}
\value{
List with components
  \item{clustID}{vector of assigned clusters}
  \item{clustres}{Result of hclust}
  \item{method}{Method inputed}
  \item{metric}{Metric used}
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
# use iris data
HClust(log(iris[,-5]), metric="pearsonCorrelation",clabel=iris[,5])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
