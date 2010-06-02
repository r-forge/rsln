\name{sample.regime.ssu}
\alias{sample.regime.ssu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sample.regime.ssu(nregimes, mu.vector, tau.vector, pi.matrix, dat, current.x, ndat = length(dat))
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
  \item{pi.matrix}{
%%     ~~Describe \code{pi.matrix} here~~
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
function(nregimes,mu.vector,tau.vector,pi.matrix,dat,current.x,ndat=length(dat)){
out <- numeric(ndat)
out[1] <- sample.int(nregimes,1,prob=dnorm(dat[1],mu.vector,1/sqrt(tau.vector))*pi.matrix[,current.x[2]]*get.stationary.dist(pi.matrix))
for(i in 2:(ndat-1)) out[i] <- sample.int(nregimes,1,prob=dnorm(dat[i],mu.vector,1/sqrt(tau.vector))*pi.matrix[out[i-1],]*pi.matrix[,current.x[i+1]])
out[ndat] <- sample.int(nregimes,1,prob=dnorm(dat[ndat],mu.vector,1/sqrt(tau.vector))*pi.matrix[out[ndat-1],])
out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
