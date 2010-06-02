sample.regime.ssuI <-
function(nregimes,pi.matrix,dat,x.vec,prior.mu,prior.n,prior.alpha,prior.beta,ndat=length(dat)){
out <- x.vec
out[1] <- sample.int(nregimes,1,prob=dfun.no.mp(dat,1,out,nregimes,prior.mu,prior.n,prior.alpha,prior.beta)*pi.matrix[,x.vec[2]]*get.stationary.dist(pi.matrix))
for(i in 2:(ndat-1)) out[i] <- sample.int(nregimes,1,prob=dfun.no.mp(dat,i,out,nregimes,prior.mu,prior.n,prior.alpha,prior.beta)*pi.matrix[out[i-1],]*pi.matrix[,x.vec[i+1]])
out[ndat] <- sample.int(nregimes,1,prob=dfun.no.mp(dat,ndat,out,nregimes,prior.mu,prior.n,prior.alpha,prior.beta)*pi.matrix[out[ndat-1],])
out
}

