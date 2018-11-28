\name{criterion.matrix}
\alias{criterion.matrix}
\docType{data}
\title{
Function to calculate the cutoff criterion given vectors of ratios and z-score or p-values
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{
criterion.matrix(zscore,rat,pval,hirat=1.2,lowrat=0.83,zlim=1.5,
                            plim=0.05,USEZSCORE=T)
}
\format{
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\arguments{
	\item{zscore}{
	
	}
	\item{rat}{
%%     ~~Describe \code{rat} here~~
}
  \item{pval}{
%%     ~~Describe \code{pval} here~~
}
  
  \item{hirat}{
%%     ~~Describe \code{hirat} here~~
}
  \item{lowrat}{
%%     ~~Describe \code{lowrat} here~~
}
\item{zlim}{
%%     ~~Describe \code{plim} here~~
}
  \item{plim}{
%%     ~~Describe \code{plim} here~~
}
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{

## The function currently defined as

criterion.matrix = function(zscore,rat,pval,hirat=1.2,lowrat=0.83,zlim=1.5,
                            plim=0.05,USEZSCORE=T)
{
  if(USEZSCORE) {
    crit.res = matrix(NA, nrow(zscore), ncol(zscore))
    
    notNA = !is.na(zscore) & !is.na(rat)
    
    crit.res[notNA] = 0
    
    crit.res[notNA & rat > hirat & abs(zscore) > zlim] = 1 ## up regulated
    
    crit.res[notNA & rat < lowrat & abs(zscore) > zlim] = -1 ## down regulated
  } else { ## pval
    crit.res = matrix(NA, nrow(pval), ncol(pval))
    
    notNA = !is.na(pval) & !is.na(rat)
    
    crit.res[notNA] = 0
    
    crit.res[notNA & rat > hirat & pval < plim] = 1 ## up regulated
    
    crit.res[notNA & rat < lowrat & pval < plim] = -1 ## down regulated
  }
  crit.res
}
## maybe str(criterion.matrix) ; plot(criterion.matrix) ...
}
\keyword{datasets}
