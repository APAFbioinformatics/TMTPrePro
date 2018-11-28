

PairedComp = function(dat, Label1,Label2,...)
{
    if(!(length(Label1)) == length(Label2)){
      stop("The length of two labels are not the same!")
    }

    ex.res = extractLabels(dat, Type="Paired", Label1,Label2)

    dat.res = ex.res$dat

  results.list = list()
  
  keep = max(5,grep("\\w+\\.w+$",colnames(dat.res)))

  dat.keep = dat.res[,c(1:keep)]
  
  rat.list = list()
  
  for (i in 1:length(Label1)){
    
    cols = grep(paste(Label1[[i]],".",Label2[[i]],sep=""),
                colnames(dat.res))    
    
    if(!(length(cols) == 3)) stop(paste(Label1[[i]],"to",
                                        Label2[[i]],"ratio, count or varability could NOT
                                  be extracted"))   

      
    dat.cols = dat.res[,cols]
    
    
    rat = dat.cols[,grep(paste(Label1[[i]],".",Label2[[i]],"$",sep=""),
                         names(dat.cols))]
    count = dat.cols[,grep(paste(Label1[[i]],".",Label2[[i]],".Count",sep=""),
                           names(dat.cols))]
    var = dat.cols[,grep(paste(Label1[[i]],".",Label2[[i]],".Variability",sep=""),
                         names(dat.cols))]
    
    zscore = get.zscore(rat, var, count) 
    pval <- zscore
    
    crit = criterion.matrix(zscore,rat,...,USEZSCORE=T)
    
    results.list[[i]] = data.frame(dat.cols,Zscore=zscore,Class=crit)
    rat.list[[i]] = data.frame(rat)
  }


  rr=Reduce(function(x,y){data.frame(x,y)}, rat.list)
  
  LogMeans = apply(log(rr),1,FUN=function(v) {mean(na.omit(v))})
  
  Means = exp(LogMeans)  ## Geometric mean
  
  results = data.frame(dat.keep,results.list,LogMeans,Means)

  data.frame(results)
}




