\name{printStyle}
\alias{printStyle}
\title{
Function to print to an excel spreadsheet using styles provided by XLConnect
}
\description{
}
\usage{
printStyle(dat, ratios, pvals, wb, tabName="results",  hiCutoff = 1.5, lowCutoff=0.67, pvalCutoff=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
	data frame to print
}
  \item{ratios}{
	Indices of ratio columns to be highlighted as up and down regulated
}
  \item{pvals}{
	Indices of ratio columns to be highlighted as up and down regulated
}
  \item{wb}{
	Workbook object to print to
}
  \item{tabName}{
	Excel tab name to generate
}
   \item{hiCutoff}{
	High threshold for highlighting ratios in pink
}
   \item{lowCutoff}{
	Low threshold for highlighting ratios in green
}
   \item{pvalCutoff}{
	Threshold for highlighting pvalues in yellow
}

}
\details{

}
\value{
Spreadsheet with highlighting printed
}
\references{
}
\author{
D. Pascovici
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link{XLConnect}}
}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
