\name{extractLabels}
\alias{extractLabels}
\docType{data}
\title{
Function to extract specific labels from an TMT proteome discoverer file
}
\description{

}
\usage{
extractLabels(dat, Label1, Label2, keep = min(10, ncol(dat)))	
}
\arguments{
  \item{dat}{
	data frame containing TMT ratios and others
}
  \item{Label1}{
	Collection of TMT labels; these will be the numerators: Label1.Label2
}
  \item{Label2}{
 	Collection of TMT labels; these will be the denominators: Label1.Label2
}
  \item{keep}{
	Number of starting columns to keep, in addition to the mapped labels
}
}
\value{
List with 2 items:
	\item{dat}{
	data frame with columns extracted correspoding to labels from Label1.Label2
	}
	\item{N}{
	Total number of columns found
	}
}
\format{
  The format is:
 chr "extractLabels"
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
dat = read.dpiscoverer(fname)

ex.res = extractLabels(dat, c("127_N 127_C"),c("126 128_C"))

}
\keyword{datasets}
