get.zscore=
  function(rat,var,count)
  {   
    df = data.frame(cbind(rat, var, count))
    
    names(df) = c("rat","var","cnt")   
    
    df$zscore = NA
    
    if(length(which(df$var==0))>0){
      df[which(df$var==0),"var"] = 0.001 ## adjust the 0
    }
    
    df[!is.na(df$var),"zscore"] = 
      100*log(df[!is.na(df$var),"rat"])/df[!is.na(df$var),"var"]
    
    dat.res = data.frame(zscore=df$zscore)
     
    dat.res
  }