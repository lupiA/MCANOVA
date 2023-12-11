#' Finds disjoint chromosome segments of a minimum length and size
#'
#' @description Calls getSegmentsChr function to compute small chromosome segments
#' 
#' @param bp ...
#' base pair position
#' @param chr ...
#' chromosome
#' @param minBPSize ...
#' minimum base pair size for each segment
#' @param minSize ...
#' minimum number of SNPs for each segment
#' @param firstSegment ...
#' default 1
#' @param verbose ...
#' default FALSE
#'
#' @return A numeric vector of segment labels corresponding to each marker
#'
#' @export

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
