rsln.fit <-
function(data.values,nregimes,burnin=0,ndraws=1000,prior.mu=0,prior.n=0.001,prior.alpha=0.001,prior.beta=0.001,
prior.dirichlet=matrix(1,nregimes,nregimes),pi.start=matrix(1/nregimes,nregimes,nregimes),state.sampler=c("ssu","ffbs")[1],constrained=FALSE){
dat <- log(data.values)
#Set up recording
ndat <- length(dat)
pi.mat <- matrix(NA,ndraws,nregimes^2)
current.pi <- pi.start
mu.mat <- tau.mat <- matrix(0,ndraws,nregimes)
x.mat <- matrix(NA,ndraws,ndat)
x.mat[1,] <- sample.int(nregimes,ndat,replace=T,prob=get.stationary.dist(pi.start))
while(mean(order(tapply(dat,x.mat[1,],mean))==1:nregimes)<1) x.mat[1,] <- sample.int(nregimes,ndat,replace=T,prob=get.stationary.dist(pi.start))
lag.count <- get.lag.count(x.mat[1,],nregimes)
for(i in 1:nregimes) current.pi[i,] <- rdirichlet(prior.dirichlet[i,] + lag.count[i,])
while(mean(order(mu.mat[1,])==1:nregimes)<1|mu.mat[1,1]==0){
reg.param.out <- regime.param(dat,x.mat[1,],nregimes,prior.mu,prior.n,prior.alpha,prior.beta)
mu.mat[1,] <- reg.param.out$mu
tau.mat[1,] <- reg.param.out$tau
}
# BURN IN ###########################################
for(rp in 1:burnin){
lag.count <- get.lag.count(x.mat[1,],nregimes)
for(i in 1:nregimes) current.pi[i,] <- rdirichlet(prior.dirichlet[i,] + lag.count[i,])
reg.param.out <- regime.param(dat,x.mat[1,],nregimes,prior.mu,prior.n,prior.alpha,prior.beta)
if(constrained&mean(order(reg.param.out$mu)==1:nregimes)<1) {
} else {
mu.mat[1,] <- reg.param.out$mu
tau.mat[1,] <- reg.param.out$tau
}
if(state.sampler=="ssu") x.mat[1,] <- sample.regime.ssu(nregimes,mu.mat[1,],tau.mat[1,],current.pi,dat,x.mat[1,],ndat)
if(state.sampler=="ffbs") x.mat[1,] <- sample.regime.ffbs(nregimes,mu.mat[1,],tau.mat[1,],current.pi,dat,x.mat[1,],ndat)
if(state.sampler=="ssuI") x.mat[1,] <- sample.regime.ssuI(nregimes,current.pi,dat,x.mat[1,],prior.mu,prior.n,prior.alpha,prior.beta)
if(state.sampler=="ffbsI") x.mat[1,] <- sample.regime.ffbsI(nregimes,current.pi,dat,x.mat[1,],prior.mu,prior.n,prior.alpha,prior.beta)
}
# SAMPLER ###########################################
for(rp in 1:ndraws){
lag.count <- get.lag.count(x.mat[rp,],nregimes)
for(i in 1:nregimes) current.pi[i,] <- rdirichlet(prior.dirichlet[i,] + lag.count[i,])
pi.mat[rp,] <- as.vector(current.pi)
reg.param.out <- regime.param(dat,x.mat[rp,],nregimes,prior.mu,prior.n,prior.alpha,prior.beta)
if(constrained&mean(order(reg.param.out$mu)==1:nregimes)<1) {
mu.mat[rp,] <- mu.mat[max(1,rp-1),]
tau.mat[rp,] <- tau.mat[max(1,rp-1),]
} else {
mu.mat[rp,] <- reg.param.out$mu
tau.mat[rp,] <- reg.param.out$tau
}
if(state.sampler=="ssu") if(rp!=ndraws) x.mat[rp+1,] <- sample.regime.ssu(nregimes,mu.mat[rp,],tau.mat[rp,],current.pi,dat,x.mat[rp,])
if(state.sampler=="ffbs") if(rp!=ndraws) x.mat[rp+1,] <- sample.regime.ffbs(nregimes,mu.mat[rp,],tau.mat[rp,],current.pi,dat,x.mat[rp,])
if(state.sampler=="ssuI") if(rp!=ndraws) x.mat[rp+1,] <- sample.regime.ssuI(nregimes,current.pi,dat,x.mat[rp,],prior.mu,prior.n,prior.alpha,prior.beta)
if(state.sampler=="ffbsI") if(rp!=ndraws) x.mat[rp+1,] <- sample.regime.ffbsI(nregimes,current.pi,dat,x.mat[rp,],prior.mu,prior.n,prior.alpha,prior.beta)
}
list(regime.mat=x.mat,mu.mat=mu.mat,tau.mat=tau.mat,sig2=1/tau.mat,pi.mat=pi.mat)
}

