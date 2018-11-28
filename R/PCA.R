PCA <-
function (data, labelValue, scaleR = FALSE, scaleC = TRUE, k = min(dim(data)) - 
    1) 
{
    if (k > min(dim(data) - 1)) 
        warning("The number of components was too large compared to the data and was adjusted accordingly")
    k <- min(k, min(dim(data)) - 1)
    if (scaleR) {
        row.nrm <- apply(data, 1, sd)
        row.nrm <- pmax(row.nrm, 1e-04)
        data <- sweep(data, 1, row.nrm, FUN = "/")
    }
    result <- try(prcomp(data, retx = TRUE, scale = scaleC), 
        silent = TRUE)
    if (inherits(result, "try-error")) 
        stop("Failed to Calculate Principal Components")
    componentVariances <- result$sdev^2
    componentLoadings <- result$rotation[, 1:k]
    componentScores <- result$x[, 1:k]
    totalVariance <- sum(componentVariances)
    componentVariances <- componentVariances[1:k]
    z <- componentScores
    plot(cloud(z[, 1] ~ z[, 3] + z[, 2], groups = as.factor(labelValue), 
        auto.key = list(points = TRUE, pch = 19, space = "right"), 
        xlab = "PC 3", ylab = "PC 2", zlab = "PC 1", distance = 0.1, 
        main = "Projection in the space of the first 3 princial components"))
    value <- list(componentVariances = componentVariances, componentScores = componentScores, 
        componentLoadings = componentLoadings, summary = summary(result))
    value
}
