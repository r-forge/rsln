rsln.gen <-
function(ndata,nregimes,mean.vec=rep(0,nregimes),prec.vec=rep(1,nregimes),trans.mat=matrix(1/nregimes,nregimes,nregimes)){
#Error Checking
if(length(mean.vec)!=nregimes) stop("Length of mean vector must equal number of regimes")
if(length(prec.vec)!=nregimes) stop("Length of precision vector must equal number of regimes")
if(!is.matrix(trans.mat)) stop("Transition matrix must be a matrix")
if(sum(dim(trans.mat)!=nregimes)>0) stop("Dimension of transition matrix must be RxR (R is the number of regimes)")
reg.values <- numeric(ndata)
reg.values[1] <- sample.int(nregimes,1,prob=get.stationary.dist(trans.mat))
for(i in 2:ndata) reg.values[i] <- sample.int(nregimes,1,prob=trans.mat[reg.values[i-1],])
data.values <- exp(sapply(reg.values,function(xx) rnorm(1,mean.vec[xx],1/prec.vec[xx])))
list(data.values=data.values,regimes=reg.values,true.mean=mean.vec,true.prec=prec.vec,true.pi=trans.mat)
}

