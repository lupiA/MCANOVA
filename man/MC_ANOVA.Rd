\name{MC_ANOVA}
\alias{MC_ANOVA}
\title{Performs MC-ANOVA}
\description{
This function predicts genetic values drawn from the core using SNPs not in the core and those in the core that are not randomly chosen to be QTL
}
\usage{
MC_ANOVA(X, X2 = NULL, core, nQTL, nRep = NULL, maxRep = 300, lambda = 1e-8, sampler = rnorm, ...)
}
\arguments{
    \item{X}{genotype matrix for group 1}
    \item{X2}{genotype matrix for group 2}
    \item{core}{columns of X that represent the core for which RSq is estimated}
    \item{nQTL}{number of causal loci for the MC simulation}
    \item{nRep}{number of monte carlo replicates, if null, it is internally determined}
    \item{maxRep}{maximum number of MC replicates}
    \item{lambda}{shrinkage parameter}
    \item{sampler}{the function to be used to sample effects}
}
\value{
    A matrix with ancestry groups in the rows, MC-ANOVA correlation estimates in the first column, and Monte Carlo error estimates in the second column
}
