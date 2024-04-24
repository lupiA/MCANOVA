GD.R <- function(XX,Xy,p=ncol(XX),b=rep(0,p), nIter=50,learning_rate=rep(.1,p),lambda=0,b0=rep(0,p),lambda0=1,returnPath=FALSE,map_flag=FALSE){

    learning_rate=learning_rate/mean(diag(XX))
    LR = rep(NA,p)
    B=matrix(nrow=ncol(XX),ncol=nIter,NA)
    B[,1]=b
    for(i in 2:ncol(B)){
      b=B[,i-1]
      for(j in 1:p){
        gradient=sum(XX[,j]*b)+lambda*b[j]-Xy[j]-lambda*lambda0*b0[j]
        b[j]=B[j,i-1]-learning_rate[j]*gradient
        LR[j]=learning_rate[j]*gradient
      }
      B[,i]=b
      print(i)
    }
    rownames(B)=rownames(XX)
    colnames(B)=paste0('iter_',1:(nIter))
    if(!returnPath){
      B=as.vector(B[,ncol(B)])
    }
    return(list(B,LR,learning_rate))
}
