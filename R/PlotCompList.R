



PlotCompList = function(res.list,comparisonInfo,...){
  
  for(jj in 1:length(res.list))
  {
    if(names(res.list)[[jj]] != "AllData"){
      
      comp = names(res.list)[[jj]]
      
      ratios=grep("1\\w+$",colnames(res.list[[jj]]))
      
      ratname = paste(comparisonInfo[jj,3],comparisonInfo[jj,4],sep="/")
      

      zscores=grep("^zscore", colnames(res.list[[jj]]))
      
      means=grep("^Means$", colnames(res.list[[jj]]))

      ttestpval=grep("^oneSplTTest$",colnames(res.list[[jj]]))
      
      if((length(means))==1 & (length(ttestpval)) ==1){ # one sample t-test
        png(paste(comp,"1TTestVolcano.png"),2000,2000,res=300)
        plotVolcano(res.list[[jj]][,ttestpval], res.list[[jj]][,means], 
                    labels=res.list[[jj]][,"Accession"],ratname,...,USEZSCORE=F)
        dev.off()
      } else if (length(ratios)==1 & length(zscores)==1){ # single
        
        png(paste(comp, "Volcano.png"), 2000, 2000, res=300)
        plotVolcano(res.list[[jj]][,zscores], res.list[[jj]][,ratios], 
                    labels=res.list[[jj]][,"Accession"],ratname,...)
        dev.off()  
      } else { # Paired
        crit.col <- grep("^Class", names(res.list[[jj]]))
        crit.mat <- res.list[[jj]][,crit.col, drop=FALSE]
        crit.idx <- apply( (!is.na(crit.mat)) & (abs(crit.mat) > 0), 1, FUN=sum)
        common.idx <- which(crit.idx> 0)
        
        if (sum(common.idx > 0)) {
          png(paste(comp, "Paired.png", sep=""), 1000, 1000, res=150)
          plot(log(res.list[[jj]][common.idx,ratios]),
               main=paste(ratname,"Paired Comparison",sep="\n"))
          abline(v=0)
          abline(h=0)
          dev.off()  
      } 
    } # End paired
  }
}
}

