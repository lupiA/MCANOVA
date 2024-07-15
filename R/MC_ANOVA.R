#' Performs MC-ANOVA
#'
#' @description This function predicts genetic values drawn from the core using SNPs not in the core
#' and those in the core that are not randomly chosen to be QTL
#'
#' @param X 
#' genotype matrix for group 1
#' @param X2 ...
#' genotype matrix for group 2
#' @param core ...
#' columns of X that represent the core for which RSq is estimated
#' @param nQTL ...
#' number of causal loci for the MC simulation
#' @param nRep ...
#' number of monte carlo replicates, if null, it is internally determined
#' @param maxRep ...
#' maximum number of MC replicates
#' @param lambda ...
#' shrinkage parameter
#' @param sampler ...
#' the function to be used to sample effects
#'
#' @return A matrix with ancestry groups in the rows
#'
#' @export

MC_ANOVA <- function(X, X2 = NULL, core, nQTL, nRep = NULL, maxRep = 300, lambda = 1e-8, sampler = rnorm, weights = NULL, ...) {

    pop2 <- !is.null(X2)
    if (pop2) {
        stopifnot(ncol(X) == ncol(X2))
    }

    if (is.character(core)) {
        core <- colnames(X) %in% core
    }
    if (is.logical(core)) {
        core <- which(core)
    }

    X <- preprocess(X, impute = TRUE, center = TRUE, scale = FALSE)
    C <- crossprod(X)
    diag(C) <- diag(C) + lambda

    # Compute CINV
    if (pop2) {
        X2 <- preprocess(X2, impute = TRUE, center = TRUE, scale = FALSE)
    }

    nCore <- length(core)
    if (nCore < nQTL) { nQTL <- nCore }

    if (nCore == 1) {
        b <- solve(C[-core, -core, drop = FALSE], C[-core, core])

        COR <- cor(X[, core], X[, -core, drop = FALSE] %*% b)
        if (pop2) {
            gHat2 <- X2[, -core, drop = FALSE] %*% b
            COR <- c('pop1' = COR, pop2 = cor(X2[, core], gHat2))
        }
    } else {
        if (is.null(nRep)) {
            nRep <- min(maxRep, choose(nCore, nQTL) * 100)
        }
        COR <- matrix(nrow = nRep, ncol = ifelse(pop2, 2, 1), NA)

        for (i in 1:nRep) {
            a <- sampler(n = nQTL, ...)

            if(!is.null(weights)){
              cols <- sort(sample(core, size = nQTL, prob = weights, replace = FALSE))
            } else{
              cols <- sort(sample(core, size = nQTL))
            }
            
            rhs <- C[-cols, cols] %*% a
            C11 <- C[-cols, -cols, drop = FALSE]
            C12 <- C[-cols, cols, drop = FALSE]
            C22 <- C[cols, cols, drop = FALSE]
            INV <- chol2inv(chol(C11))
            b <- crossprod(INV, rhs)

            SS0 <- t(a) %*% C22 %*% a
            COND.VAR <- C22 - crossprod(C12, crossprod(INV, C12))
            SS1 <- t(a) %*% COND.VAR %*% a
            RSq <- 1 - SS1 / SS0

            COR[i, 1] <- sqrt(RSq)
            if (pop2) {
                g2 <- X2[, cols, drop = FALSE] %*% a
                gHat2 <- X2[, -cols, drop = FALSE] %*% b

                COR[i, 2] <- cor(g2, gHat2)
            }
        }
    }

    ANS=matrix(nrow = ifelse(pop2, 2, 1),ncol = 6, NA)
    ANS[, 1] <- colMeans(COR)
    ANS[, 2] <- apply(FUN = sd, X = COR, MARGIN = 2)
    ANS[, 3] <- apply(FUN = sd, X = COR^2, MARGIN = 2)
    ANS[, 4] <- mean(RA.res, na.rm = T)
    ANS[, 5] <- apply(FUN = sd, X = RA.res, MARGIN = 2)
    ANS[, 6] <- nrow(COR)
    colnames(ANS) <- c('Cor', 'Cor_SD', 'Rsq_SD', 'RA','RA_SD','nRep')

    if (pop2) {
        rownames(ANS) <- c('Group 1', 'Group 2')
    } else {
        rownames(ANS) <- c('Group 1')
    }
    return(ANS)
}
