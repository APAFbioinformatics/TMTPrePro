\name{read.pdiscoverer}
\alias{read.pdiscoverer}
\docType{data}
\title{
Read a proteome discoverer protein group file
}
\description{
}
\usage{
read.pdiscoverer(fname,cutoof)
}
\arguments{
  \item{fname}{
	Proteome discoverer protein group file name
}
  \item{cutoff}{
	TO DEFINE
}
 
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\value{
 Data frame containing protein records
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
# sample file from package
fname <- paste(system.file("files", package="TMTRpac"),"SetABC_RedeinedDB_proteingroups.txt", sep="/")
dat = read.pdiscoverer(pdFile)

}
\keyword{datasets}
