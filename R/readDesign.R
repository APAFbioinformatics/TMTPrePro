


######################################################################################
# file format: Label, Replicate (optional), Group (optional), Comparion type, ....
######################################################################################




readDesign <- function(designFile)
{
  
  design <- readWorkbook(designFile,1)  #Label, Replicates, Group, comp type x, ... 

  # check the row names
  
  if ( !tolower("Label") %in% tolower(names(design)) )
    stop("The design file must contain Label column!")
  
  
  comparisonInfo <- data.frame( Name=character(), Type=character(), Numerator.Label=character(), Denominator.Label=character(),
                     stringsAsFactors=F)
  
  
  
  
  # Extract single comparisons
  
   
    idx <- grep("single", tolower(names(design)) )
	
	if(length(idx) > 0) {
    
    for(ii in 1:length(idx) ) {
      
      top <- paste(design[grep("t", tolower(design[,idx[ii]])), 1], collapse = " ")
      
      bot <- paste(design[grep("b", tolower(design[,idx[ii]])), 1], collapse = " ")
      
      comparisonInfo <- rbind(comparisonInfo, 
                              data.frame(Name=paste("SingleComp", ii), 
                                         Type="Single", Numerator.Label=top, Denominator.Label=bot) )
      
    }
    
    
  }
  
   
  
  # Extract onesampe t test comparisons
  
    
    idx <- grep("onesample", tolower(names(design)) )
	
	if(length(idx) > 0) {
    
    for(ii in 1:length(idx) ) {
      
      top <- paste(design[grep("t", tolower(design[,idx[ii]])), 1], collapse = " ")
                   
      bot <- paste(design[grep("b", tolower(design[,idx[ii]])), 1], collapse = " ")
      
      comparisonInfo <- rbind(comparisonInfo, 
                   data.frame(Name=paste("OneSampleTTestComp", ii), 
                              Type="OneSampleTTest", Numerator.Label=top, Denominator.Label=bot) )
    
    }
  }
  
  
  
  # Extract paired comparisons
  
  
    idx <- grep("pair", tolower(names(design)) )
    
	if(length(idx) > 0) {
	
    for(ii in 1:length(idx) ) { # (t1, t2) (b1, b2)
      
      top1 <- paste(design[grep("t1", tolower(design[,idx[ii]])), 1], collapse = " ")
      top2 <- paste(design[grep("t2", tolower(design[,idx[ii]])), 1], collapse = " ")
      
      bot1 <- paste(design[grep("b1", tolower(design[,idx[ii]])), 1], collapse = " ")
      bot2 <- paste(design[grep("b2", tolower(design[,idx[ii]])), 1], collapse = " ")
      
  
      
      if ( !all.equal(sapply(list(top1,top2,bot1,bot2), length), rep(1,4)) )          
        stop("T1, T2, B1, B2 must only have one occurrances in Paired comparison")
      
      comparisonInfo <- rbind(comparisonInfo, 
                   data.frame(Name=paste("PairedComp", ii), 
                              Type="Paired", Numerator.Label=paste(top1,top2, sep=" "), Denominator.Label=paste(bot1,bot2, sep=" ")) )
      
    }
  }
  
 
  
  comparisonInfo  #Name, Type, Numerator.Label, Denominator.Label

  
}







