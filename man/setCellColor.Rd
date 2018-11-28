\name{setCellColor}
\alias{setCellColor}
\docType{data}
\title{
Set an xlsx cell color
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{
setCellColor(wb,sheet,row.num,col.num,col="yellow")
}
\arguments{
wb - a workbook
sheet - a worksheet
row.num - the row number of the cell
col.num - the column number of the cell
col - the color to be set
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
library(xlsx)

wb = createWorkbook()
sheet = createSheet(wb, "test")
##fname <- paste(system.file("files", package="TMTRpac"),"SetABC_RedeinedDB_proteingroups.txt", sep="/")

##tmt <- read.pdiscoverer(fname)
						
##addDataFrame(tmt,sheet,row.names=F)

row.num=6
col.num=9
setCellColor(wb,sheet,row.num,col.num,"lightblue")

saveWorkbook(wb,"test.xlsx")
}
\keyword{datasets}
