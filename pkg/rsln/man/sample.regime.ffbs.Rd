\name{sample.regime.ffbs}
\alias{sample.regime.ffbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sample.regime.ffbs(nregimes, mu.vector, tau.vector, current.pi, dat, current.x, ndat = length(dat))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nregimes}{
%%     ~~Describe \code{nregimes} here~~
}
  \item{mu.vector}{
%%     ~~Describe \code{mu.vector} here~~
}
  \item{tau.vector}{
%%     ~~Describe \code{tau.vector} here~~
}
  \item{current.pi}{
%%     ~~Describe \code{current.pi} here~~
}
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{current.x}{
%%     ~~Describe \code{current.x} here~~
}
  \item{ndat}{
%%     ~~Describe \code{ndat} here~~
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
