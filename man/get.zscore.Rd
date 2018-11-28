\name{get.zscore}
\alias{get.zscore}
\docType{data}
\title{
Get the Z score from the log ratio and varability
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{
get.zscore(ratio,var)
}
\format{
ratio - a vector of the protein ratio
var - a vector of the varability
count - a vector of the count
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
ratios = c(0.732,0.874,1.027)
vars = c(56.3,21.0,22.1)
counts = c(12,30,63)
z.score=get.zscore(ratios,vars,counts) 
	
}
\keyword{datasets}
