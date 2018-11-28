
extractLabels <- function(dat, Type, Label1, Label2, keep=min(6, ncol(dat)) )  {  
  
  Label1 <- unique(Label1)
  Label2 <- unique(Label2)  
  
  to.ext = switch(Type, 
					Single=as.vector(paste(Label1, Label2)),
					Paired=as.vector(paste(Label1, Label2)),
					OneSampleTTest=as.vector(outer(Label1,  Label2, FUN=paste)) )
  
  to.ext.1 <- gsub("(.*) (.*)", "\\1", to.ext)
  to.ext.2 <- gsub("(.*) (.*)", "\\2", to.ext)
  
  idx <- (to.ext.1 != to.ext.2)
  
  to.ext <- to.ext[idx]
  to.ext <- gsub(" ", "\\.", to.ext)    
  
  lab.exist<- unlist(lapply(to.ext, FUN=function(x)
    { length(grep(tolower(x), tolower(names(dat))))}))
  
  lab.notextidx<- which(lab.exist==0)
  
  if(length(lab.notextidx) != 0)
    warnings(paste(paste(to.ext[lab.notextidx],collapse =","), 
               "do NOT exist! 
               Please check the labels in the comparison file.",sep=" "))
  
  lab.notexist <- data.frame(to.ext[lab.notextidx])
  
  cols <- unlist(lapply(to.ext, FUN=function(x)
    { grep(tolower(x), tolower(names(dat))) }))
  
  
  # keep at most up to "Peptides"
  
  keep <- min(keep, grep("Peptides\\.", names(dat)))
  
  
  list(dat=dat[, c(1:keep, cols)], N=length(cols))  
  
}