


SingleComp = function(dat,Label1,Label2,...)
  
{
  ## extract columns
    
  ex.res = extractLabels(dat,Type="Single", Label1,Label2)
 
  dat.res=ex.res$dat
  
  
  if(ex.res$N==3){
    
    var = dat.res[,grep("variability",tolower(names(dat.res))),drop=FALSE]
 
    rat = dat.res[,grep(tolower(paste(Label1,".",Label2,
                                      "$",sep="")),tolower(names(dat.res))),drop=FALSE]
    
    count = dat.res[,grep("count",tolower(names(dat.res))),drop=FALSE]
    
    zscore = get.zscore(rat, var, count) 
    pval <- zscore
    
    crit = criterion.matrix(zscore,rat,pval,...,USEZSCORE=T)
  
    dat.res = data.frame(dat.res,Zscore=zscore,Class=crit)
  } else {
    stop("The requested labels for the single comparision were not found in the dataset")
    
  }
  
  dat.res
}