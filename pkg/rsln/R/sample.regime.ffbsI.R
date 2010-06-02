sample.regime.ffbsI <-
function(nregimes,current.pi,dat,x.vec,prior.mu,prior.n,prior.alpha,prior.beta,ndat=length(dat)){
ck <- numeric(ndat)
phi.mat <- matrix(NA,ndat,nregimes)
phi.temp <- get.stationary.dist(current.pi)
beta.mat <- matrix(NA,ndat,nregimes)
out <- x.vec
for(i in 1:ndat){
ck[i] <- sum(phi.temp*dfun.no.mp(dat,i,out,nregimes,prior.mu,prior.n,prior.alpha,prior.beta))
phi.mat[i,] <- phi.temp*dfun.no.mp(dat,i,out,nregimes,prior.mu,prior.n,prior.alpha,prior.beta)/ck[i]
phi.temp <- apply(phi.mat[i,]*t(current.pi),2,sum)
}
out[ndat] <- sample.int(nregimes,1,prob=phi.mat[ndat,])
for(i in ((ndat-1):1)){
temp <- phi.mat[i,]*current.pi[,out[i+1]]
out[i] <- sample.int(nregimes,1,prob=temp)
}
out
}

