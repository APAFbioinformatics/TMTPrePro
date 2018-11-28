apply.t1 = function(v, log=T)
{
  res = NA
  
  v=na.omit(v)
  if(log) v=log(v)
  
  if(length(v) > 1){
    tp = try(t.test(v)$p.value)
    if(!inherits(tp, "try-error")) res = tp
  }
  
  res
}