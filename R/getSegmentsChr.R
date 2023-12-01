getSegmentsChr <- function(bp, minBPSize = 10e3, minSize = 10, firstSegment = 1) {
    x <- as.character(cut(bp, breaks = seq(from = 0, to = max(bp, na.rm = TRUE) + minBPSize + 10, by = minBPSize)))
    y <- as.integer(factor(x, levels = unique(x)))
    segmentsY <- unique(y)
    nSegments <- length(segmentsY)

    p <- length(y)
    segment <- rep(0, p)
    currentSegment <- 1
    nextSegment <- firstSegment

    while (currentSegment <= nSegments) {
        tmp <- which(y == segmentsY[currentSegment])
        tmpLength <- length(tmp)

        while ((tmpLength < minSize) && (currentSegment < nSegments)) {
            currentSegment <- currentSegment + 1
            tmp <- c(tmp, which(y == segmentsY[currentSegment]))
            tmpLength <- length(tmp)
        }

        segment[tmp] <- nextSegment
        currentSegment <- currentSegment + 1
        nextSegment <- nextSegment + 1
    }

    tmp <- sum(segment == max(segment, na.rm = TRUE))
    if (tmp < minSize) {
        segment[segment == max(segment, na.rm = TRUE)] <- max(segment, na.rm = TRUE) - 1
    }

    return(segment)
}
