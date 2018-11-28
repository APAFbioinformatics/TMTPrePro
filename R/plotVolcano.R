plotVolcano <-
  function(pzval, rat, labels, ratname,hirat=1.2,lowrat=0.83,zlim=1.5,plim=0.05,
           USEZSCORE=T){    
    if(USEZSCORE) {# use zscore
      
      print.idx <- ( (log(abs(pzval)) > log(zlim)) & 
                       ( (rat < lowrat) | (rat > hirat) ) )
      
      plot(t(log(rat)), t(log(abs(pzval))), pch=".", cex=1.5,
           xlab="Ln ratio", ylab="Ln(abs(Z Score))",
           main=paste(ratname,"Ratio VS Z-score","Single Comparison",sep="\n"))
      abline(h = log(zlim), col =  "red")
      abline(v =  log(hirat),    col = "blue")
      abline(v = log(lowrat),    col = "blue")      
        
      text(t(log(rat[print.idx])), log(abs(pzval))[print.idx], cex=0.3, col="red",
           labels=labels[print.idx])
      
    } else { # use pvalue
      
    print.idx <- ( pzval < plim & ( (rat < lowrat) | (rat > hirat) ) )
      
    plot(t(log(rat)), -log(t(pzval)), pch=".", cex=1.5, 
         xlab="Ln ratio", ylab="- Ln p value",
         main=paste(ratname,"Ratio VS p value","One Sample T-test",sep="\n"))
    abline(h = -log(plim), col =  "red")
    abline(v =  log(hirat),    col = "blue")
    abline(v = log(lowrat),    col = "blue")    
       
    text(t(log(rat[print.idx])), -log(pzval)[print.idx], cex=0.3, col="red",
         labels=labels[print.idx])    
    
  }
}