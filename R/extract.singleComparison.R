## extract values for a single comparison
extract.singleComparison=
  function(tmt.df,numera_label,denom_label)
  {
    label.name=paste(numera_label,denom_label,sep=".")
    col.nums=grep(label.name,colnames(tmt.df))
    ## 3 columns: ratio, count, variability
    
    start.col=min(grep("1",colnames(tmt.df)))
    extract.cols=c(1:(start.col-1))
    extract.cols=c(extract.cols,col.nums) 
    
    tmt.df[,extract.cols]
  }