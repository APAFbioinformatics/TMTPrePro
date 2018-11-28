`plotDensity` <- function (data, group=rownames(data), xlab="Abundance") 
{
group <- as.factor(group)
colours <- rainbow(length(levels(group)))
col <- colours[group]

    # par(mfrow = c(1, 1))
    x <- t(as.matrix(data))
    ndx <- rep(FALSE, ncol(x))
    for (i in 1:ncol(x)) ndx[i] <- is.numeric(x[, i])
    if (sum(ndx) > 0) {
        x <- x[, ndx, drop = FALSE]
        try(d1 <- density(as.vector(as.matrix(x))))
        if (inherits(d1, "try-error")) 
            Error("Failed to generate the density plots")
        ymx <- max(d1$y)
        plot(d1, type = "n", xlab = xlab , ylab = "Density", xlim=c(-2,2),
            main = "Sample Densities", ylim = c(0, 2 * ymx), 
            yaxp = c(0, 2 * ymx, 5))
        for (i in 1:ncol(x)) {
            try(d1 <- density(x[, i]))
            if (inherits(d1, "try-error")) 
                Error(paste("Failed to generate the density plot for sample", 
                  i))
            lines(d1, lty = i, col = col[i])
        }
 legend("topright", legend=levels(group), col=colours, lty=1)
    }
}
