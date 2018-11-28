


OneSampleTComp <- function(dat,Label1,Label2,...)
  
  {    
    ## extract columns
    
    ex.res = extractLabels(dat, Type="OneSampleTTest", Label1,Label2)  
    
    dat.res=ex.res$dat
 
    cols <- grep(tolower("1\\w+$"),tolower(colnames(dat.res))) ## ratio
    
    rat=dat.res[,cols, drop=F]
    
     
    LogMeans = apply(log(rat), 1, FUN=function(v) {mean(na.omit(v))})
    Means = exp(LogMeans)
    
    ## One sample t test for log ratios
    
    oneSplTTest = apply(rat,1,FUN=apply.t1)
    TTestAdjusted = p.adjust(oneSplTTest, method="fdr" )  

    
    crit = criterion.matrix(zscore=as.data.frame(oneSplTTest),  ## dummy
                            rat=as.data.frame(Means),
                            pval=as.data.frame(oneSplTTest),
                            ...,
                            USEZSCORE=F)
    
    dat.res=data.frame(dat.res, Means, LogMeans, oneSplTTest, Class=crit)
    
    dat.res[order(dat.res$oneSplTTest),]
  }


