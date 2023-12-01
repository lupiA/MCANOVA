getSegments <- function(bp, chr, minBPSize = 10e3, minSize = 10, firstSegment = 1, verbose = FALSE) {
    p <- length(bp)
    segments <- rep(NA, p)
    
    chromosomes <- unique(chr)
    nextSegment <- 1
    
    for (i in 1:length(chromosomes)) {
        if (verbose) {
            message('Chr ', chromosomes[i])
        }
        tmp <- which(chr == chromosomes[i])
        
        chr_wind <- getSegmentsChr(bp[tmp], minBPSize = minBPSize, minSize = minSize, firstSegment = nextSegment)
        segments[tmp] <- chr_wind
        nextSegment <- max(segments, na.rm = TRUE) + 1
    }
    return(segments)
}
