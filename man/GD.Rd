\name{GD}
\alias{GD}
\title{Runs a gradient descent algorithm}
\description{
This function performs a gradient descent with the option to have SNP-specific learning rates.
}
\usage{
GD(XX,Xy,p=ncol(XX),b=rep(0,p),nIter=50,learning_rate=rep(.1,p),lambda=0,b0=rep(0,p),lambda0=1,returnPath=FALSE,map_flag=FALSE)
}
\arguments{
    \item{XX}{lhs}
    \item{Xy}{rhs}
    \item{p}{number of columns for lhs matrix}
    \item{b}{initial values (e.g., EU PGS)}
    \item{nIter}{number of iterations}
    \item{learning_rate}{learning rate; can be fixed or adaptive (supply as vector)}
    \item{lambda}{shrinkage parameter}
    \item{b0}{additional shrinkage parameter}
    \item{lambda0}{value to shrink towards}
    \item{returnPath}{true or false}
    \item{map_flag}{true or false}
}
\value{
    List of estimates, the learning rate multiplied by gradient, and the learning rate.
}
