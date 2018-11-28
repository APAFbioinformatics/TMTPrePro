updateComparisonInfo <- function(comparisonInfo, dat.err)
{

  comparisonInfo <- comparisonInfo[!dat.err$comp.notexist,]
  
  comparisonInfo
}