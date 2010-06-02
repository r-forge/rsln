sample.regime.ssu <-
function(nregimes,mu.vector,tau.vector,pi.matrix,dat,current.x,ndat=length(dat)){
out <- numeric(ndat)
out[1] <- sample.int(nregimes,1,prob=dnorm(dat[1],mu.vector,1/sqrt(tau.vector))*pi.matrix[,current.x[2]]*get.stationary.dist(pi.matrix))
for(i in 2:(ndat-1)) out[i] <- sample.int(nregimes,1,prob=dnorm(dat[i],mu.vector,1/sqrt(tau.vector))*pi.matrix[out[i-1],]*pi.matrix[,current.x[i+1]])
out[ndat] <- sample.int(nregimes,1,prob=dnorm(dat[ndat],mu.vector,1/sqrt(tau.vector))*pi.matrix[out[ndat-1],])
out
}

