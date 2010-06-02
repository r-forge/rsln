regime.param <-
function(dat,current.x,nregimes,mu.0,gamma.0,alpha.0,beta.0){
mu.out <- tau.out <- numeric(nregimes)
for(reg in 1:nregimes){
regime.data <- dat[current.x==reg]
if(length(regime.data)==0) regime.data <- dat
regime.length <- length(regime.data)
mu.n <- (mu.0*gamma.0 + regime.length * mean(regime.data))/(gamma.0 + regime.length)
gamma.n <- gamma.0 + regime.length
alpha.n <- alpha.0 + regime.length/2
beta.n <- beta.0 + sum((regime.data-mean(regime.data))^2)/2 + (regime.length*gamma.0*(mean(regime.data)-mu.0)^2)/(2*(regime.length + gamma.0))
tau.out[reg] <- rgamma(1,alpha.n,beta.n)
mu.out[reg] <- rnorm(1,mu.n, sqrt(1/(gamma.n*tau.out[reg])))
}
list(mu=mu.out,tau=tau.out)
}

