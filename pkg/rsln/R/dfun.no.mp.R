dfun.no.mp <-
function(dat,index,x.vec,nregimes,prior.mu,prior.n,prior.alpha,prior.beta){
datmi <- dat[-index]; xmi <- x.vec[-index]
x.bar <- tapply(datmi,xmi,mean)
s2 <- numeric(nregimes)
for(i in 1:nregimes) s2[i] <- mean((datmi[xmi==i] - x.bar[i])^2)
nn <- tapply(datmi,xmi,length)
mu.n <- (prior.n + nn)^(-1)*(prior.n*prior.mu + nn*x.bar)
beta.n <- prior.beta + nn*s2/2 + prior.n*nn*(prior.mu-x.bar)^2/(2*(nn+prior.n))
tdensity(dat[index],mu.n,(nn+prior.n)/(nn + prior.n + 1)*(prior.alpha + nn/2)/beta.n,2*prior.alpha + nn)
}

