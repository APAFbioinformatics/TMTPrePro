cordist <- function(x) { as.dist((1-cor(t(x)))/2) }
