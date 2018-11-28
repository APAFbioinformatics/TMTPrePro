


printCompList = function(res.list, comparisonInfo, dat.err,
                         labelNames=NULL, excelName="Results.xlsx")
{
  wb = createWorkbook()
  
  
  
  addWorksheet(wb,"Comparison")
  
  
  writeData(wb, sheet = "Comparison", x = comparisonInfo, rowNames = F)
  
  
  addWorksheet(wb, "ErrorNote")
  
  
  writeData(wb, sheet = "ErrorNote", x = dat.err$lab.notexist, rowNames = F)
  
  ## print each single run comparison 
 
  upStyle <- createStyle(fgFill = "yellow")
  downStyle <- createStyle(fgFill = "#AFB9BA")
  
  for(i in 1:length(res.list))    
  {
    res.list[[i]][is.na(res.list[[i]])] <- NA
    
    Sheet <- names(res.list)[[i]]
    addWorksheet(wb, names(res.list)[[i]]) 

    writeData(wb, sheet = Sheet, x= res.list[[i]], rowNames = F)
    
    if(names(res.list)[[i]] != "AllData") { # only highlight for comparisions 
      
      crit.cols = grep("Class",colnames(res.list[[i]]))
      
       
      if("Zscore" %in% tolower(names(res.list[[i]]))) { # single  and paired    
      
        rat.col = grep(tolower("1\\w+$"),tolower(names(res.list[[i]])))
        
        score.col = grep(tolower("Zscore"),tolower(names(res.list[[i]])))
        
        if(!(length(rat.col) == length(crit.cols))) stop("Number of ratios and class are different")
         
        for(ii in 1:length(crit.cols)){
          
          up.rows <- which(res.list[[i]][,crit.cols[ii]] == 1)
          down.rows <- which(res.list[[i]][,crit.cols[ii]] == -1)
          
          addStyle(wb, Sheet, style = upStyle, rows = up.rows+1, cols = rat.col, gridExpand = T)
          addStyle(wb, Sheet, style = upStyle, rows = up.rows+1, cols = score.col, gridExpand = T)
          
          addStyle(wb, Sheet, style = downStyle, rows = down.rows+1, cols = rat.col, gridExpand = T)
          addStyle(wb, Sheet, style = downStyle, rows = down.rows+1, cols = score.col, gridExpand = T)
      
        }          
      
      } else if ("onesplttest" %in% tolower(names(res.list[[i]]))){ # one sample t-test 
         
        rat.col = grep(tolower("^means"),tolower(names(res.list[[i]])))
        
        score.col = grep(tolower("ttest$"),tolower(names(res.list[[i]]))) 
        
        if(!(length(rat.col) == length(crit.cols))) stop("Number of ratios and crit are different")
        
        for(ii in 1:length(crit.cols)){
          up.rows <- which(res.list[[i]][,crit.cols[ii]] == 1)
          down.rows <- which(res.list[[i]][,crit.cols[ii]] == -1)
          
          addStyle(wb, Sheet, style = upStyle, rows = up.rows+1, cols = rat.col, gridExpand = T)
          addStyle(wb, Sheet, style = upStyle, rows = up.rows+1, cols = score.col, gridExpand = T)
          
          addStyle(wb, Sheet, style = downStyle, rows = down.rows+1, cols = rat.col, gridExpand = T)
          addStyle(wb, Sheet, style = downStyle, rows = down.rows+1, cols = score.col, gridExpand = T)
          
        }  
      }
            
    }
  }
  
  saveWorkbook(wb,excelName, overwrite = T)
}


