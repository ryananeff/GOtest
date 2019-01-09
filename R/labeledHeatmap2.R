#modified from labeledHeatmap in package WGCNA
labeledHeatmap2=function (Matrix, xLabels, yLabels = NULL, xSymbols = NULL, ySymbols = NULL, 
    colorLabels = NULL, xColorLabels = FALSE, yColorLabels = FALSE, 
    checkColorsValid = TRUE, invertColors = FALSE, setStdMargins = TRUE, 
    xLabelsPosition = "bottom", xLabelsAngle = 45, xLabelsAdj = 1, 
    xColorWidth = 0.05, yColorWidth = 0.05, xColorOffset = par("cxy")[1]/3, 
    yColorOffset = par("cxy")[2]/3, colors = NULL, naColor = "grey", 
    textMatrix = NULL, cex.text = NULL, textAdj = c(0.5, 0.5), 
    cex.lab = NULL, cex.lab.x = cex.lab, cex.lab.y = cex.lab, cex.legend=cex.lab,
    colors.lab.x = 1, colors.lab.y = 1, bg.lab.x = NULL, bg.lab.y = NULL, 
    plotLegend = TRUE, keepLegendSpace = plotLegend, verticalSeparator.x = NULL, 
    verticalSeparator.col = 1, verticalSeparator.lty = 1, verticalSeparator.lwd = 1, 
    verticalSeparator.ext = 0, horizontalSeparator.y = NULL, 
    horizontalSeparator.col = 1, horizontalSeparator.lty = 1, 
    horizontalSeparator.lwd = 1, horizontalSeparator.ext = 0, 
    ...) 
{
    if (!is.null(colorLabels)) {
        xColorLabels = colorLabels
        yColorLabels = colorLabels
    }
    if (is.null(yLabels) & (!is.null(xLabels)) & (dim(Matrix)[1] ==  dim(Matrix)[2])) yLabels = xLabels
    nCols = ncol(Matrix)
    nRows = nrow(Matrix)
    if (checkColorsValid) {
        xValidColors = !is.na(match(substring(xLabels, 3), colors()))
        yValidColors = !is.na(match(substring(yLabels, 3), colors()))
    }
    else {
        xValidColors = rep(TRUE, length(xLabels))
        yValidColors = rep(TRUE, length(yLabels))
    }
    if (sum(xValidColors) > 0) 
        xColorLabInd = c(1:length(xLabels))[xValidColors]
    if (sum(!xValidColors) > 0) 
        xTextLabInd = c(1:length(xLabels))[!xValidColors]
    if (sum(yValidColors) > 0) 
        yColorLabInd = c(1:length(yLabels))[yValidColors]
    if (sum(!yValidColors) > 0) 
        yTextLabInd = c(1:length(yLabels))[!yValidColors]
    if (setStdMargins) {
        if (xColorLabels & yColorLabels) {
            par(mar = c(2, 2, 3, 5) + 0.2)
        }
        else {
            par(mar = c(7, 7, 3, 5) + 0.2)
        }
    }
    xLabPos = charmatch(xLabelsPosition, c("bottom", "top"))
    if (is.na(xLabPos)) stop("Argument 'xLabelsPosition' must be (a unique abbreviation of) 'bottom', 'top'")
    if (is.null(colors)) colors = heat.colors(30)
    if (invertColors) colors = rev(colors)
    labPos = heatmapWithLegend2(Matrix, signed = FALSE, colors = colors, 
        naColor = naColor, cex.legend = cex.legend, plotLegend = plotLegend, 
        keepLegendSpace = keepLegendSpace, ...)
    nxlabels = length(xLabels)
    plotbox = labPos$box
    xmin = plotbox[1]
    xmax = plotbox[2]
    ymin = plotbox[3]
    yrange = plotbox[4] - ymin
    ymax = plotbox[4]
    xrange = xmax - xmin
    xLeft = labPos$xLeft
    xRight = labPos$xRight
    yTop = labPos$yTop
    yBot = labPos$yBot
    xspacing = labPos$xMid[2] - labPos$xMid[1]
    yspacing = abs(labPos$yMid[2] - labPos$yMid[1])
    nylabels = length(yLabels)
    offsetx = yColorOffset
    offsety = xColorOffset
    xColW = min(xmax - xmin, ymax - ymin) * xColorWidth
    yColW = min(xmax - xmin, ymax - ymin) * yColorWidth
    if (any(xValidColors)) 
        offsety = offsety + xColW
    if (any(yValidColors)) 
        offsetx = offsetx + yColW
    extension.left = par("mai")[2] * par("cxy")[1]/par("cin")[1]
    extension.bottom = par("mai")[1] * par("cxy")[2]/par("cin")[2] - 
        offsety
    extension.top = par("mai")[3] * par("cxy")[2]/par("cin")[2] - 
        offsety
    figureBox = par("usr")
    figXrange = figureBox[2] - figureBox[1]
    figYrange = figureBox[4] - figureBox[3]
    if (!is.null(bg.lab.x)) {
        bg.lab.x = extend2(bg.lab.x, nCols)
        if (xLabPos == 1) {
            y0 = ymin
            ext = extension.bottom
            sign = 1
        }
        else {
            y0 = ymax
            ext = extension.top
            sign = -1
        }
        figureDims = par("pin")
        angle = xLabelsAngle/180 * pi
        ratio = figureDims[1]/figureDims[2] * figYrange/figXrange
        ext.x = -sign * ext * 1/tan(angle)/ratio
        ext.y = sign * ext * sign(sin(angle))
        for (c in 1:nCols) polygon(x = c(xLeft[c], xLeft[c], 
            xLeft[c] + ext.x, xRight[c] + ext.x, xRight[c], xRight[c]), 
            y = c(y0, y0 - sign * offsety, y0 - sign * offsety - 
                ext.y, y0 - sign * offsety - ext.y, y0 - sign * 
                offsety, y0), border = bg.lab.x[c], col = bg.lab.x[c], 
            xpd = TRUE)
    }
    if (!is.null(bg.lab.y)) {
        bg.lab.y = extend2(bg.lab.y, nRows)
        reverseRows = TRUE
        if (reverseRows) {
            bg.lab.y = rev(bg.lab.y)
        }
        for (r in 1:nRows) rect(xmin - extension.left, yBot[r], 
            xmin, yTop[r], col = bg.lab.y[r], border = bg.lab.y[r], 
            xpd = TRUE)
    }
    if (sum(!xValidColors) > 0) {
        xLabYPos = ifelse(xLabPos == 1, ymin - offsety, ymax + 
            offsety)
        if (is.null(cex.lab)) 
            cex.lab = 1
        text(labPos$xMid[xTextLabInd], xLabYPos, srt = xLabelsAngle, 
            adj = xLabelsAdj, labels = xLabels[xTextLabInd], 
            xpd = TRUE, cex = cex.lab.x, col = colors.lab.x)
    }
    if (sum(xValidColors) > 0) {
        baseY = ifelse(xLabPos == 1, ymin - offsety, ymax + offsety)
        deltaY = ifelse(xLabPos == 1, xColW, -xColW)
        rect(xleft = labPos$xMid[xColorLabInd] - xspacing/2, 
            ybottom = baseY, xright = labPos$xMid[xColorLabInd] + 
                xspacing/2, ytop = baseY + deltaY, density = -1, 
            col = substring(xLabels[xColorLabInd], 3), border = substring(xLabels[xColorLabInd], 
                3), xpd = TRUE)
        if (!is.null(xSymbols)) 
            text(labPos$xMid[xColorLabInd], baseY - sign(deltaY) * 
                offsety, xSymbols[xColorLabInd], adj = xLabelsAdj, 
                xpd = TRUE, srt = xLabelsAngle, cex = cex.lab.x, 
                col = colors.lab.x)
    }
    if (sum(!yValidColors) > 0) {
        if (is.null(cex.lab)) 
            cex.lab = 1
        text(xmin - offsetx, labPos$yMid[yTextLabInd], srt = 0, 
            adj = c(1, 0.5), labels = yLabels[yTextLabInd], xpd = TRUE, 
            cex = cex.lab.y, col = colors.lab.y)
    }
    if (sum(yValidColors) > 0) {
        rect(xleft = xmin - offsetx, ybottom = rev(labPos$yMid[yColorLabInd]) - 
            yspacing/2, xright = xmin - offsetx + yColW, ytop = rev(labPos$yMid[yColorLabInd]) + 
            yspacing/2, density = -1, col = substring(rev(yLabels[yColorLabInd]), 
            3), border = substring(rev(yLabels[yColorLabInd]), 
            3), xpd = TRUE)
        if (!is.null(ySymbols)) 
            text(xmin + yColW - 2 * offsetx, labPos$yMid[yColorLabInd], 
                ySymbols[yColorLabInd], adj = c(1, 0.5), xpd = TRUE, 
                cex = cex.lab.y, col = colors.lab.y)
    }
    if (length(verticalSeparator.x) > 0) {
        nLines = length(verticalSeparator.x)
        vs.col = extend2(verticalSeparator.col, nLines)
        vs.lty = extend2(verticalSeparator.lty, nLines)
        vs.lwd = extend2(verticalSeparator.lwd, nLines)
        vs.ext = extend2(verticalSeparator.ext, nLines)
        if (any(verticalSeparator.x < 0 | verticalSeparator.x > 
            nCols)) 
            stop("If given. 'verticalSeparator.x' must all be between 0 and the number of columns.")
        x.lines = ifelse(verticalSeparator.x > 0, labPos$xRight[verticalSeparator.x], 
            labPos$xLeft[1])
        for (l in 1:nLines) lines(rep(x.lines[l], 2), c(ymin, 
            ymax), col = vs.col[l], lty = vs.lty[l], lwd = vs.lwd[l])
        angle = xLabelsAngle/180 * pi
        if (xLabelsPosition == "bottom") {
            sign = 1
            y0 = ymin
        }
        else {
            sign = -1
            y0 = ymax
        }
        figureDims = par("pin")
        ratio = figureDims[1]/figureDims[2] * figYrange/figXrange
        ext.x = -sign * extension.bottom * 1/tan(angle)/ratio
        ext.y = sign * extension.bottom * sign(sin(angle))
        for (l in 1:nLines) lines(c(x.lines[l], x.lines[l], x.lines[l] + 
            vs.ext * ext.x), c(y0, y0 - sign * offsety, y0 - 
            sign * offsety - vs.ext * ext.y), col = vs.col[l], 
            lty = vs.lty[l], lwd = vs.lwd[l], xpd = TRUE)
    }
    if (length(horizontalSeparator.y) > 0) {
        if (any(horizontalSeparator.y < 0 | horizontalSeparator.y > 
            nRows)) 
            stop("If given. 'horizontalSeparator.y' must all be between 0 and the number of rows.")
        reverseRows = TRUE
        if (reverseRows) {
            horizontalSeparator.y = nRows - horizontalSeparator.y + 
                1
            y.lines = ifelse(horizontalSeparator.y <= nRows, 
                labPos$yBot[horizontalSeparator.y], labPos$yTop[nRows])
        }
        else {
            y.lines = ifelse(horizontalSeparator.y > 0, labPos$yBot[horizontalSeparator.y], 
                labPos$yTop[1])
        }
        nLines = length(horizontalSeparator.y)
        vs.col = extend2(horizontalSeparator.col, nLines)
        vs.lty = extend2(horizontalSeparator.lty, nLines)
        vs.lwd = extend2(horizontalSeparator.lwd, nLines)
        vs.ext = extend2(horizontalSeparator.ext, nLines)
        for (l in 1:nLines) lines(c(xmin - vs.ext[l] * extension.left, 
            xmax), rep(y.lines[l], 2), col = vs.col[l], lty = vs.lty[l], 
            lwd = vs.lwd[l], xpd = TRUE)
    }
    if (!is.null(textMatrix)) {
        if (is.null(cex.text)) 
            cex.text = par("cex")
        if (is.null(dim(textMatrix))) 
            if (length(textMatrix) == prod(dim(Matrix))) 
                dim(textMatrix) = dim(Matrix)
        if (!isTRUE(all.equal(dim(textMatrix), dim(Matrix)))) 
            stop("labeledHeatmap: textMatrix was given, but has dimensions incompatible with Matrix.")
        for (rw in 1:dim(Matrix)[1]) for (cl in 1:dim(Matrix)[2]) {
            text(labPos$xMid[cl], labPos$yMid[rw], as.character(textMatrix[rw, 
                cl]), xpd = TRUE, cex = cex.text, adj = textAdj)
        }
    }
    axis(1, labels = FALSE, tick = FALSE)
    axis(2, labels = FALSE, tick = FALSE)
    axis(3, labels = FALSE, tick = FALSE)
    axis(4, labels = FALSE, tick = FALSE)
    invisible(labPos)
}
#adapted package WGCNA
blueWhiteRed2=function (n, gamma = 1, endSaturation = 1) 
{
    if (endSaturation > 1 | endSaturation < 0) 
        stop("'endSaturation' must be between 0 and 1.")
    es = 1 - endSaturation
    blueEnd = c(0.05 + es * 0.45, 0.55 + es * 0.25, 1)
    redEnd = c(1, 0.2 + es * 0.6, 0.6 * es)
    middle = c(1, 1, 1)
    half = as.integer(n/2)
    if (n%%2 == 0) {
        index1 = c(1:half)
        index2 = c(1:half) + half
        frac1 = ((index1 - 1)/(half - 1))^(1/gamma)
        frac2 = rev(frac1)
    }
    else {
        index1 = c(1:(half + 1))
        index2 = c(1:half) + half + 1
        frac1 = (c(0:half)/half)^(1/gamma)
        frac2 = rev((c(1:half)/half)^(1/gamma))
    }
    cols = matrix(0, n, 3)
    for (c in 1:3) {
        cols[index1, c] = blueEnd[c] + (middle[c] - blueEnd[c]) * 
            frac1
        cols[index2, c] = redEnd[c] + (middle[c] - redEnd[c]) * 
            frac2
    }
    rgb(cols[, 1], cols[, 2], cols[, 3], maxColorValue = 1)
}
extend2=function (x, n){
    nRep = ceiling(n/length(x))
    rep(x, nRep)[1:n]
}
heatmapWithLegend2=function (data, signed, colors, naColor = "grey", zlim = NULL, 
    reverseRows = TRUE, plotLegend = TRUE, keepLegendSpace = plotLegend, 
    cex.legend = 1, legendShrink = 0.94, legendSpace = 0.1, legendWidth = 0.02, 
    legendGap = 0.02, frame = TRUE, frameTicks = FALSE, tickLen = 0.02, 
    ...) 
{
    data = as.matrix(data)
    nCols = ncol(data)
    nRows = nrow(data)
    if (is.null(zlim)) {
        zlim = range(data, na.rm = TRUE)
        if (signed) zlim = c(-max(abs(zlim)), max(abs(zlim)))
    }
    barplot(1, col = "white", border = "white", axisnames = FALSE, axes = FALSE, ...)
    box = par("usr")
    xminAll = box[1]
    xmaxAll = box[2]
    yminAll = box[3]
    ymaxAll = box[4]
    if (!keepLegendSpace && !plotLegend) {
        legendSpace = 0
        legendWidth = 0
        legendGap = 0
    }
    ymin = yminAll
    ymax = ymaxAll
    xmin = xminAll
    xmax = xmaxAll - legendSpace * (xmaxAll - xminAll)
    xStep = (xmax - xmin)/nCols
    xLeft = xmin + c(0:(nCols - 1)) * xStep
    xRight = xLeft + xStep
    xMid = (xLeft + xRight)/2
    yStep = (ymax - ymin)/nRows
    yBot = ymin + c(0:(nRows - 1)) * yStep
    yTop = yBot + yStep
    yMid = c(yTop + yBot)/2
    if (reverseRows) {
        colorMat = numbers2colors(reverseRows2(data), signed, 
            colors = colors, lim = zlim, naColor = naColor)
    }
    else colorMat = numbers2colors(data, signed, colors = colors, 
        lim = zlim, naColor = naColor)
    dim(colorMat) = dim(data)
#	colorMat=as.raster(colorMat)
#	rasterImage(colorMat, xmin, yBot, xmax, yTop)
#	browser()
    for (c in 1:nCols) {
        rect(xleft = rep(xLeft[c], nRows), xright = rep(xRight[c], 
            nRows), ybottom = yBot, ytop = yTop, col = colorMat[, c], border = colorMat[, c])
    }
    if (frame) 
        lines(c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin))
    if (plotLegend) {
        yminL = yminAll + (1 - legendShrink) * (ymaxAll - yminAll)
        ymaxL = ymaxAll - (1 - legendShrink) * (ymaxAll - yminAll)
        xminL = xmaxAll - (xmaxAll - xminAll) * (legendSpace - 
            legendGap)
        xmaxL = xmaxAll - (xmaxAll - xminAll) * (legendSpace - 
            legendGap - legendWidth)
        tickVal = autoTicks2(zlim[1], zlim[2])
        tickY = (tickVal - zlim[1])/(zlim[2] - zlim[1]) * (ymaxL - 
            yminL) + yminL
        nTicks = length(tickVal)
        for (t in 1:nTicks) lines(c(xmaxL, xmaxL + (xmaxAll - 
            xminAll) * 0.8 * tickLen), c(tickY[t], tickY[t]))
        text(rep(xmaxL + (xmaxAll - xminAll) * tickLen), tickY, 
            tickVal, adj = c(0, 0.5), cex = cex.legend, xpd = TRUE)
        nColors = length(colors)
        ybl = (ymaxL - yminL)/nColors * (0:(nColors - 1)) + yminL
        ytl = (ymaxL - yminL)/nColors * (1:nColors) + yminL
        rect(xleft = rep(xminL, nColors), xright = rep(xmaxL, 
            nColors), ybottom = ybl, ytop = ytl, col = colors, 
            border = colors)
        lines(c(xminL, xmaxL, xmaxL, xminL, xminL), c(yminL, 
            yminL, ymaxL, ymaxL, yminL))
    }
    list(xMid = xMid, yMid = if (reverseRows) rev(yMid) else yMid, 
        box = c(xmin, xmax, ymin, ymax), xLeft = xLeft, xRight = xRight, 
        yTop = yTop, yBot = yBot)
}
reverseRows2=function (Matrix) 
{
    ind = seq(from = dim(Matrix)[1], to = 1, by = -1)
    Matrix[ind, ]
}
numbers2colors=function (x, signed = NULL, centered = signed, lim = NULL, commonLim = FALSE, 
    colors = if (signed) blueWhiteRed2(100) else blueWhiteRed2(100)[51:100], 
    naColor = "grey") 
{
    x = as.matrix(x)
    if (!is.numeric(x)) 
        stop("'x' must be numeric. For a factor, please use as.numeric(x) in the call.")
    if (is.null(signed)) {
        if (any(x < 0, na.rm = TRUE) & any(x > 0, na.rm = TRUE)) {
            signed = TRUE
        }
        else signed = FALSE
    }
    if (is.null(centered)) 
        centered = signed
    if (is.null(lim)) {
        if (signed & centered) {
            max = apply(abs(x), 2, max, na.rm = TRUE)
            lim = as.matrix(cbind(-max, max))
        }
        else {
            lim = as.matrix(cbind(apply(x, 2, min, na.rm = TRUE), 
                apply(x, 2, max, na.rm = TRUE)))
        }
        if (commonLim) 
            lim = c(min(lim[, 1], na.rm = TRUE), max(lim[, 2], 
                na.rm = TRUE))
    }
    if (is.null(dim(lim))) {
        if (length(lim) != 2) 
            stop("'lim' must be a vector of length 2 or a matrix with 2 columns.")
        if (!is.numeric(lim)) 
            stop("'lim' must be numeric")
        if (sum(is.finite(lim)) != 2) 
            stop("'lim' must be finite.")
        lim = t(as.matrix(lim))
    }
    else {
        if (ncol(x) != nrow(lim)) 
            stop("Incompatible numbers of columns in 'x' and rows in 'lim'.")
        if (!is.numeric(lim)) 
            stop("'lim' must be numeric")
        if (sum(is.finite(lim)) != length(lim)) 
            stop("'lim' must be finite.")
    }
    xMin = matrix(lim[, 1], nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
    xMax = matrix(lim[, 2], nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
    if (sum(xMin == xMax) > 0) 
        warning("(some columns in) 'x' are constant. Their color will be the color of NA.")
    xx = x
    xx[is.na(xx)] = ((xMin + xMax)[is.na(xx)])/2
    if (sum(x < xMin, na.rm = TRUE) > 0) {
        warning("Some values of 'x' are below given minimum and will be truncated to the minimum.")
        x[xx < xMin] = xMin[xx < xMin]
    }
    if (sum(x > xMax, na.rm = TRUE) > 0) {
        warning("Some values of 'x' are above given maximum and will be truncated to the maximum.")
        x[xx > xMax] = xMax[xx > xMax]
    }
    mmEq = xMin == xMax
    nColors = length(colors)
    xCol = array(naColor, dim = dim(x))
    xInd = (x - xMin)/(xMax - xMin)
    xInd[xInd == 1] = 1 - 0.5/nColors
    xCol[!mmEq] = colors[as.integer(xInd[!mmEq] * nColors) + 
        1]
    xCol[is.na(xCol)] = naColor
    xCol
}
autoTicks2=function (min, max, maxTicks = 6, tickPos = c(1, 2, 5)) 
{
    range = max - min
    tick0 = range/maxTicks
    maxTick = max(tickPos)
    mult = 1
    if (tick0 < maxTick/10) {
        while (tick0 < maxTick/10) {
            tick0 = 10 * tick0
            mult = mult * 10
        }
    }
    else while (tick0 >= maxTick) {
        tick0 = tick0/10
        mult = mult/10
    }
    ind = sum(tick0 > tickPos) + 1
    tickStep = tickPos[ind]/mult
    lowTick = min/tickStep
    if (floor(lowTick) != lowTick) lowTick = lowTick + 1
    lowTick = floor(lowTick)
    ticks = tickStep * (lowTick:(lowTick + maxTicks + 1))
    ticks = ticks[ticks <= max]
    ticks
}
