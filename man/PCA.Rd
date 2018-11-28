\name{PCA}
\alias{PCA}
\title{
 Function to do principal components analysis; prcomp wrapper
}
\description{

}
\usage{
PCA(data, labelValue, scaleR = FALSE, scaleC = TRUE, k = min(dim(data)) - 1)
}
\arguments{
  \item{data}{
	Data matrix, variables as columns, samples as rows
}
  \item{labelValue}{
	Factor containing the sample groups; length = nrow(data)
}
  \item{scaleR}{
	TRUE/FALSE: Should the rows (samples) be scaled? Default FALSE
}
  \item{scaleC}{
	TRUE/FALSE: Should the columns (variables) be scaled? Default TRUE
}
  \item{k}{
	Number of components; if more than the rank will be adjusted to dim(data)-1
}
}
\details{
}
\value{
List with components
\item{componentVariances}{Variances; i.e. corresponding eigenvalues squared}
\item{componentScores}{sample scores for each component}
\item{componentLoadings}{variable loadings for each component}
\item{summary}{prcomp summary, including percentage of variance explained}
}
\references{
}
\author{
}
\note{
}
\seealso{
\code{\link{prcomp}}
}
\examples{

library(lattice)

PCA(iris[,-5], iris[,5])

}

