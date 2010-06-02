sample.regime.ffbs <-
function(nregimes,mu.vector,tau.vector,current.pi,dat,current.x,ndat=length(dat)){
ck <- numeric(ndat)
phi.mat <- matrix(NA,ndat,nregimes)
phi.temp <- get.stationary.dist(current.pi)
beta.mat <- matrix(NA,ndat,nregimes)
out <- numeric(ndat)
for(i in 1:ndat){
ck[i] <- sum(phi.temp*dnorm(dat[i],mu.vector,1/sqrt(tau.vector)))
phi.mat[i,] <- phi.temp*dnorm(dat[i],mu.vector,1/sqrt(tau.vector))/ck[i]
phi.temp <- apply(phi.mat[i,]*t(current.pi),2,sum)
}
out[ndat] <- sample.int(nregimes,1,prob=phi.mat[ndat,])
for(i in ((ndat-1):1)){
temp <- phi.mat[i,]*current.pi[,out[i+1]]
out[i] <- sample.int(nregimes,1,prob=temp)
}
out
}

