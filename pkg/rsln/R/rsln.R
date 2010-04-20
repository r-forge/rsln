rm(list=ls(all=TRUE))

rsln.gen <- function(ndata,nregimes,mean.vec=rep(0,nregimes),prec.vec=rep(1,nregimes),trans.mat=matrix(1/nregimes,nregimes,nregimes)){
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
get.lag.count <- function(current.x,nregimes) as.matrix(table(factor(current.x[-length(current.x)],levels=1:nregimes),factor(current.x[-1],levels=1:nregimes)))
rdirichlet <- function(alpha.vec) {
	gamma.vec <- sapply(alpha.vec, function(xx) rgamma(1,xx))
	gamma.vec/sum(gamma.vec)
}
get.stationary.dist <- function(mat) {
	e.mat <- eigen(t(mat))
	e.vec <- e.mat$vectors[,zapsmall(e.mat$values)==1]
	e.vec/sum(e.vec)
}
regime.param <- function(xx,mu.0,gamma.0,alpha.0,beta.0){
	regime.length <- length(xx)
	mu.n <- (mu.0*gamma.0 + regime.length * mean(xx))/(gamma.0 + regime.length)
	gamma.n <- gamma.0 + regime.length
	alpha.n <- alpha.0 + regime.length/2
	beta.n <- beta.0 + sum((xx-mean(xx))^2)/2 + (regime.length*gamma.0*(mean(xx)-mu.0)^2)/(2*(regime.length + gamma.0))

	tau.c <- rgamma(1,alpha.n,beta.n)
	mu.c <- rnorm(1,mu.n, sqrt(1/(gamma.n*tau.c)))
	c(mu.c,tau.c)
}
sample.regime.ssu <- function(nregimes,mu.vector,tau.vector,pi.matrix,dat,current.x,ndat=length(dat)){
	out <- numeric(ndat)
	out[1] <- sample.int(nregimes,1,prob=dnorm(dat[1],mu.vector,1/sqrt(tau.vector))*pi.matrix[,current.x[2]])
	for(i in 2:(ndat-1)) out[i] <- sample.int(nregimes,1,prob=dnorm(dat[i],mu.vector,1/sqrt(tau.vector))*pi.matrix[out[i-1],]*pi.c[,current.x[i+1]])
	out[ndat] <- sample.int(nregimes,1,prob=dnorm(dat[ndat],mu.vector,1/sqrt(tau.vector))*pi.matrix[out[ndat-1],])
	out
}
sample.regime.ffbs <- function(nregimes,mu.vector,tau.vector,current.pi,dat,current.x,ndat=length(dat)){
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
	#beta.mat[ndat,] <- 1/ck[ndat]
	#for(i in ((ndat-1):1)){
	#	beta.mat[i,] <- (1/ck[i])*apply(pi.matrix*dnorm(dat[i+1],mu.vector,1/sqrt(tau.vector))*beta.mat[i+1,],2,sum)
	#	out[i] <- sample.int(nregimes,1,prob=phi.mat[i,]*beta.mat[i,])
	#}
	#out[ndat] <- sample.int(nregimes,1,prob=pi.matrix[out[ndat-1],]*dnorm(dat[ndat],mu.vector,1/sqrt(tau.vector))*beta.mat[ndat,])
	out
}
rsln.fit <- function(data.values,nregimes,burnin=1000,ndraws=1000,prior.mu=0,prior.n=0.001,prior.alpha=1000,prior.beta=0.001,prior.dirichlet=matrix(1,nregimes,nregimes),pi.start=matrix(1/nregimes,nregimes,nregimes)){
	dat <- log(data.values)
	#Set up recording
	ndat <- length(dat)
	mu.mat <- tau.mat <- matrix(NA,ndraws,nregimes)
	x.mat <- matrix(NA,ndraws,ndat)
	x.mat[1,] <- sample.int(nregimes,ndat,replace=T,prob=get.stationary.dist(pi.start))
	pi.mat <- matrix(NA,ndraws,nregimes^2)
	current.pi <- pi.start
	#Burn in
	for(rp in 1:burnin){
		lag.count <- get.lag.count(x.mat[1,],nregimes)
		for(i in 1:nregimes) current.pi[i,] <- rdirichlet(prior.dirichlet[i,] + lag.count[i,])
		#pi.mat[rp,] <- as.vector(pi.c)
		for(i in 1:nregimes){
			regime.data <- dat[x.mat[1,]==i]
			if(length(regime.data)==0) regime.data <- dat
			reg.param.out <- regime.param(regime.data,prior.mu,prior.n,prior.alpha,prior.beta)
			mu.mat[1,i] <- reg.param.out[1]
			tau.mat[1,i] <- reg.param.out[2]
		}
		x.mat[1,] <- sample.regime.ffbs(nregimes,mu.mat[1,],tau.mat[1,],current.pi,dat,x.mat[1,],ndat)
	}
	for(rp in 1:ndraws){
		lag.count <- get.lag.count(x.mat[rp,],nregimes)
		for(i in 1:nregimes) current.pi[i,] <- rdirichlet(prior.dirichlet[i,] + lag.count[i,])
		pi.mat[rp,] <- as.vector(current.pi)
		for(i in 1:nregimes){
			regime.data <- dat[x.mat[rp,]==i]
			if(length(regime.data)==0) regime.data <- dat
			reg.param.out <- regime.param(regime.data,prior.mu,prior.n,prior.alpha,prior.beta)
			mu.mat[rp,i] <- reg.param.out[1]
			tau.mat[rp,i] <- reg.param.out[2]
		}
		if(rp!=ndraws) x.mat[rp+1,] <- sample.regime.ffbs(nregimes,mu.mat[rp,],tau.mat[rp,],current.pi,dat,x.mat[rp,])
	}
	list(regime.mat=x.mat,mu.mat=mu.mat,tau.mat=tau.mat,sig2=1/tau.mat,pi.mat=pi.mat)
}