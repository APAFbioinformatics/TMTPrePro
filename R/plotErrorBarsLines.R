plotErrorBarsLines <-
function (v, barSizes, lines, labels = NULL, col = "blue", ylim=c(min(lines), max(lines)),
    ...) 
{
    barSizes[is.na(barSizes)] <- 0
    topBars <- v + 0.5 * barSizes
    bottomBars <- v - 0.5 * barSizes
    N <- length(v)
    if (is.null(labels)) 
        labels <- 1:N
    ylims <- c(min(bottomBars, ylim[1], min(lines)), max(topBars, 
        ylim[2], max(lines)))
    par(pch = 19, xaxt = "n")
    plot(as.numeric(labels), v, ylim = ylims, col = col, type = "b", 
        lwd = 3, ...)
    par(xaxt = "s")
    my.at <- 1:N
    axis(1, at = my.at, labels = labels)
    for (i in 1:N) {
        lines(c(i, i), c(topBars[i], bottomBars[i]))
    }
    for (i in 1:ncol(lines)) {
        lines(as.numeric(labels), lines[, i], lwd = 0.5, lty = "dotted", 
            col = "gray")
    }
}
