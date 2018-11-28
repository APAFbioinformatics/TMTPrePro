\name{apply.t1}
\alias{apply.t1}
\docType{data}
\title{
Apply a one sample t-test where possible
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{
	apply.t1(v, log)
}
\format{
  
}
\arguments{
  \item{v}{
	numerical vector
}
	\item{log}{
	logic value, representing if take logrism of the vector for T test
	}
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\value{

Returns the p-value of the test only

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
v <- c(NA, rnorm(10))
apply.t1(v,log=FALSE)
}
\keyword{datasets}
