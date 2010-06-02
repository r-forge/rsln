\name{regime.param}
\alias{regime.param}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
regime.param(dat, current.x, nregimes, mu.0, gamma.0, alpha.0, beta.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{current.x}{
%%     ~~Describe \code{current.x} here~~
}
  \item{nregimes}{
%%     ~~Describe \code{nregimes} here~~
}
  \item{mu.0}{
%%     ~~Describe \code{mu.0} here~~
}
  \item{gamma.0}{
%%     ~~Describe \code{gamma.0} here~~
}
  \item{alpha.0}{
%%     ~~Describe \code{alpha.0} here~~
}
  \item{beta.0}{
%%     ~~Describe \code{beta.0} here~~
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
