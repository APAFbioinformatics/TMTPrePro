get.upregulated=
  function(tmt.df,up.ratio=1.2,zscore=1.5)
  {
    # return a data frame
    ratio.col=min(grep("1",colnames(tmt.df)))
    count.col=ratio.col+1
    var.col=ratio.col+2

    tmt.df$z.score=get.zscore(tmt.df[,ratio.col],
                                tmt.df[,var.col]) ## pecentage
    tmt.df[,ratio.col]=as.numeric(tmt.df[,ratio.col])
    
    idx=which(tmt.df[,ratio.col]>up.ratio & 
                (abs(tmt.df$z.score)>zscore))
    tmt.df[idx,]
  }