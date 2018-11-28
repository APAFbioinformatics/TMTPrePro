get.downregulated=
  function(tmt.df,down.ratio=0.83,zscore=1.5)
  {
    # return a data frame
    ratio.col=min(grep("1",colnames(tmt.df)))
    count.col=ratio.col+1
    var.col=ratio.col+2
    
    tmt.df$z.score=get.zscore(tmt.df[,ratio.col],
                              tmt.df[,var.col]) 
    tmt.df[,ratio.col]=as.numeric(tmt.df[,ratio.col])

    idx=which(tmt.df[,ratio.col]<down.ratio & 
                abs(tmt.df$z.score)>zscore)
    tmt.df[idx,]
  }