\name{rsln.fit}
\alias{rsln.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rsln.fit(data.values, nregimes, burnin = 0, ndraws = 1000, prior.mu = 0, prior.n = 0.001, prior.alpha = 0.001, prior.beta = 0.001, prior.dirichlet = matrix(1, nregimes, nregimes), pi.start = matrix(1/nregimes, nregimes, nregimes), state.sampler = c("ssu", "ffbs")[1], constrained = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.values}{
%%     ~~Describe \code{data.values} here~~
}
  \item{nregimes}{
%%     ~~Describe \code{nregimes} here~~
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
}
  \item{ndraws}{
%%     ~~Describe \code{ndraws} here~~
}
  \item{prior.mu}{
%%     ~~Describe \code{prior.mu} here~~
}
  \item{prior.n}{
%%     ~~Describe \code{prior.n} here~~
}
  \item{prior.alpha}{
%%     ~~Describe \code{prior.alpha} here~~
}
  \item{prior.beta}{
%%     ~~Describe \code{prior.beta} here~~
}
  \item{prior.dirichlet}{
%%     ~~Describe \code{prior.dirichlet} here~~
}
  \item{pi.start}{
%%     ~~Describe \code{pi.start} here~~
}
  \item{state.sampler}{
%%     ~~Describe \code{state.sampler} here~~
}
  \item{constrained}{
%%     ~~Describe \code{constrained} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
