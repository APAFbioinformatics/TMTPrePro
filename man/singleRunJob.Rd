\name{singleRunJob}
\alias{singleRunJob}
\docType{data}
\title{
Data analysis for a single TMT run
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{
singleRunJob(fname,singleRunInfo,
                up.ratio=1.2, down.ratio=0.83,
                z.score=1.5,xlsx.name="output.xlsx")
}
\format{
 
}
\arguments{
fname	The TMT protein discoverer file name
singleRunInfo	A csv file defining all comparisons of a single run
up.ratio	the cutoff of up regulated protein expression value
up.ratio	the cutoff of down regulated protein expression value
z.score		the z score cutoff of differentially expressed proteins 
xlsx.name	the output file name
}
\value{
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
fname <- paste(system.file("files", package="TMTRpac"),"SetABC_RedeinedDB_proteingroups.txt", sep="/")
singleRunInfo<-paste(system.file("files",package="TMTRpac"),"comparisonInfo.csv",sep="/")
singleRunJob(fname,singleRunInfo)
}
\keyword{datasets}
