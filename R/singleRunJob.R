
##################################
# Default cutoffs:
# uprat=1.2; downrat=0.83; zlim=1.5
# library(TMTPrePro); library(openxlsx)
##################################


singleRunJob = function(pdFile,comparisonInfoFile, ...)
{
  dat = try(read.pdiscoverer(pdFile))
  
  if (inherits(dat,"try-error")) 
    stop ("Error reading Proteome Discoverer file!")
  
  
  comparisonInfo = try(readDesign(comparisonInfoFile) )
  
  if (inherits(comparisonInfo, "try-error"))
    stop ("Error reading design file!")
  
  
  
  dat.err = try(checkLabelExistence(dat,comparisonInfo) ) ## not existing labels
  
  if (inherits(dat.err, "try-error"))
    stop ("Error checkLabelExistence!")
  
  
  
  comparisonInfoNew = try(updateComparisonInfo(comparisonInfo, dat.err)  )## remove non-existing labels
  
  if(nrow(comparisonInfoNew) == 0)
	stop("No comparison is generated due to no existing ratio labels were found! Please check the design.")
 
  if( inherits(comparisonInfoNew, "try-error"))
    stop ("Error updateComparisonInfo!")


  res.list = try(GenCompList(dat, comparisonInfoNew, ...))
  
  if (inherits(res.list,"try-error")) 
    stop("Error in generating the list of comparision")
  
  

  tp<-try(PlotCompList(res.list,comparisonInfoNew,...))
  
  if (inherits(res.list,"try-error")) 
    stop("Error in plotting the comparison list figures")
  
  
  
  
  tp<-try(printCompList(res.list, comparisonInfoNew, dat.err))
  if (inherits(res.list,"try-error")) 
    stop("Error in printing the list of comparision")
  
  
  
  
}


