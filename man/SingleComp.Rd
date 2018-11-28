\name{SingleComp}
\alias{SingleComp}
\docType{data}
\title{
Single comparison of 2 labels of TMT data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{
SingleComp(dat,Label1,Label2)
}
\arguments{
  \item{dat}{
	A dataframe of TMT dataset
}
  \item{Label1}{
	A TMT label as the numerator label
}
  \item{Label2}{
	A TMT label as the demonimator label
}
}
\value{
Return a dataframe of single comparison
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

res.list = SingleComp(dat, c("127_C"),c("126"))
}
\keyword{datasets}
