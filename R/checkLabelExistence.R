

checkLabelExistence <- function(dat, comparisonInfo)
{
  comparisonInfo = as.matrix(comparisonInfo)
  
  comp.idx=!(comparisonInfo[,1] == "")
  comp.list=comparisonInfo[,1][comp.idx]
  comp.type=comparisonInfo[,2][comp.idx]
  
  label.list = list()
  
  for(ii in 1:length(comp.list)){
    label.list[[ii]] = list(strsplit(comparisonInfo[ii,3],split="[[:space:]]")[[1]],
                            strsplit(comparisonInfo[ii,4],split="[[:space:]]")[[1]])
    
  }
  
  # generate all comparisons
  
  lab.notexist = list()
  comp.notexist = rep(FALSE, length(comp.list) )
  
  for(j in 1:length(comp.list)){
    Type = comp.type[j]
    
    Label1=unlist(label.list[[j]][1])
    Label2=unlist(label.list[[j]][2])
    
    Label1 <- unique(Label1)
    Label2 <- unique(Label2)  
    
	
    #to.ext <- as.vector(outer(Label1,  Label2, FUN=paste)) # OneSampleTTest
    
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
    
    if(length(lab.notextidx) != 0) comp.notexist[j] = TRUE
    
    lab.notexist[[j]] <- data.frame(NotExistingLabels=to.ext[lab.notextidx])
    
  }
  
  list(lab.notexist=do.call(rbind,lab.notexist), comp.notexist=comp.notexist )
}




















