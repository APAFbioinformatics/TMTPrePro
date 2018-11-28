

GenCompList = function(dat, comparisonInfo,...)
{
  results.list=list()
  
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
  
  for(j in 1:length(comp.list)){
    Type = comp.type[j]
    
    Label1=unlist(label.list[[j]][1])
    Label2=unlist(label.list[[j]][2])
    
    
    results.list[[j]]=switch(Type, Single=SingleComp(dat,Label1,Label2,...),
                          Paired=PairedComp(dat,Label1,Label2,...),
                          OneSampleTTest=OneSampleTComp(dat,Label1,Label2,...))
  
  }
  
  names(results.list) = comp.list
  results.list[["AllData"]] = dat 
  
  ## write AllData
  write.csv(dat,file="AllData.csv")
  
  results.list
}




