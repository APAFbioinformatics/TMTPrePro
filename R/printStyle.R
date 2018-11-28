printStyle  <- function (dat, ratios, pvals, wb, tabName = "results", hiCutoff = 1.5, lowCutoff=0.67, pvalCutoff=0.05) 
{

    setStyleAction(wb, XLC$STYLE_ACTION.NONE)
    createSheet(wb, name = tabName)
    wb[tabName] <- dat
    upReg <- createCellStyle(wb)
    setFillPattern(upReg, fill = XLC$FILL.SOLID_FOREGROUND)
    setFillForegroundColor(upReg, color = XLC$COLOR.ROSE)
    downReg <- createCellStyle(wb)
    setFillPattern(downReg, fill = XLC$FILL.SOLID_FOREGROUND)
    setFillForegroundColor(downReg, color = 57)
    sigStyle <- createCellStyle(wb)
    setFillPattern(sigStyle, fill = XLC$FILL.SOLID_FOREGROUND)
    setFillForegroundColor(sigStyle, color = 43)
    for (rat in ratios) {
        up.idx <- which(!is.na(dat[, rat]) & (dat[, rat] > hiCutoff))
        if (length(up.idx) > 1) 
            setCellStyle(wb, sheet = tabName, row = 1 + up.idx, 
                col = rat, cellstyle = upReg)
        down.idx <- which(!is.na(dat[, rat]) & (dat[, rat] < 
            lowCutoff))
        if (length(down.idx) > 1) 
            setCellStyle(wb, sheet = tabName, row = 1 + down.idx, 
                col = rat, cellstyle = downReg)
    }
    for (pval in pvals) {
        sig.idx <- which(!is.na(dat[, pval]) & (dat[, pval] < 
            pvalCutoff))
        if (length(sig.idx) > 1) 
            setCellStyle(wb, sheet = tabName, row = 1 + sig.idx, 
                col = pval, cellstyle = sigStyle)
    }
}
