\name{PairedComp}
\alias{PairedComp}
\docType{data}
\title{
Paired comparison based on a number of TMT labels from the same run
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{
PairedComp(dat, Label1, Label2)
}
\arguments{
  \item{dat}{
	TMT data for that run
}
  \item{Label1}{
	Set of labels - they will be the numerators for ratios
}
  \item{Label2}{
	Set of labels - they will be the denominators for ratios
}
}
\value{
Return a list of dataframe of paired comparison
}
\format{
 
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

dat.res <- PairedComp(dat,c("127_N 127_C"),c("126 128_C"))
}
\keyword{datasets}
