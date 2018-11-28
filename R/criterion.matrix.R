


criterion.matrix = function(zscore,rat,pval,hirat = 1.5, lowrat = 0.67, zlim = 1.5,
                            plim = 0.05,USEZSCORE=T)
{
  if(USEZSCORE) {
    crit.res = matrix(NA, nrow(zscore), ncol(zscore))
    
    notNA = !is.na(zscore) & !is.na(rat)
    
    crit.res[notNA] = 0
    
    crit.res[notNA & rat > hirat & abs(zscore) > zlim] = 1 ## up regulated
    
    crit.res[notNA & rat < lowrat & abs(zscore) > zlim] = -1 ## down regulated
  } else { ## pval
    crit.res = matrix(NA, nrow(pval), ncol(pval))
    
    notNA = !is.na(pval) & !is.na(rat)
    
    crit.res[notNA] = 0
    
    crit.res[notNA & rat > hirat & pval < plim] = 1 ## up regulated
    
    crit.res[notNA & rat < lowrat & pval < plim] = -1 ## down regulated
  }
  
  crit.res
}