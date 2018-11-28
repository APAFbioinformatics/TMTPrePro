\name{oneSampleTComp}
\alias{oneSampleTComp}
\docType{data}
\title{
Function to compare TMT ratios by a one sample t-test
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{
OneSampleTComp(dat, Label1, Label2)
}
\arguments{
  \item{dat}{
	A TMT data set
}
  \item{Label1}{
A set of labels as numerator labels
}
  \item{Label2}{
A set of labels as denominator labels
}
}
\value{
Return a vector of the one sample t-test p values
}
\format{
  The format is:
 chr "oneSampleTComp"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
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

dat.res <- OneSampleTTest(dat,c("127_N 127_C"),c("126 128_C"))

}
\keyword{datasets}
