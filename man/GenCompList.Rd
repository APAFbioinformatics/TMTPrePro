\name{GenCompList}
\alias{GenCompList}
\docType{data}
\title{
Function to run various comparisons of single run TMT data
}
\description{
}
\usage{
GenCompList(dat, comparisonInfo)
}

\arguments{
  \item{dat}{
	The TMT proteome discoverer data set for one run
}
  \item{comparisonInfo}{
	The matrix of comparisons
}
}
\format{
  The format is:
 chr "GenCompList"
}
\details{
The matrix of comparisons should contain 4 columns: the analysis name,
the two sets of labels used, and the type of analysis, which can be Single,
Paired or OneSampleTTest.
}

\value{
A list with each individual comparison, and the combined data.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
fname <- paste(system.file("files", package="TMTRpac"),"SetABC_RedeinedDB_proteingroups.txt", sep="/")
dat = read.pdiscoverer(pdFile)

comparisonInfo <- data.frame(Name=c("Comp1","Comp2","Comp3"),
                             Type = c("Single","OneSampleTTest","Paired"),
                             Label1 = c("127_C","127_N 127_C","127_N 127_C"),
                             Label2 = c("126","126 128_C","126 128_C"))
							 
res.list<-GenCompList(dat,comparisonInfo)
                             
}
\keyword{datasets}
