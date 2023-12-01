\name{getSegments}
\alias{getSegments}
\title{Computes small chromosome segments}
\description{
    Calls getSegmentsChr function to compute small chromosome segments
}
\usage{
getSegments(bp, chr, minBPSize = 10e3, minSize = 10, firstSegment = 1, verbose = FALSE)
}
\arguments{
    \item{bp}{base pair position}
    \item{chr}{chromosome}
    \item{minBPSize}{minimum base pair size for each segment}
    \item{minSize}{minimum number of SNPs for each segment}
    \item{firstSegment}{default 1}
    \item{verbose}{default FALSE}
}
\value{
    A numeric vector of segment labels corresponding to each marker
}
