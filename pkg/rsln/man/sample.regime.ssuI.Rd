\name{sample.regime.ssuI}
\alias{sample.regime.ssuI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sample.regime.ssuI(nregimes, pi.matrix, dat, x.vec, prior.mu, prior.n, prior.alpha, prior.beta, ndat = length(dat))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nregimes}{
%%     ~~Describe \code{nregimes} here~~
}
  \item{pi.matrix}{
%%     ~~Describe \code{pi.matrix} here~~
}
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{x.vec}{
%%     ~~Describe \code{x.vec} here~~
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
function(nregimes,pi.matrix,dat,x.vec,prior.mu,prior.n,prior.alpha,prior.beta,ndat=length(dat)){
out <- x.vec
out[1] <- sample.int(nregimes,1,prob=dfun.no.mp(dat,1,out,nregimes,prior.mu,prior.n,prior.alpha,prior.beta)*pi.matrix[,x.vec[2]]*get.stationary.dist(pi.matrix))
for(i in 2:(ndat-1)) out[i] <- sample.int(nregimes,1,prob=dfun.no.mp(dat,i,out,nregimes,prior.mu,prior.n,prior.alpha,prior.beta)*pi.matrix[out[i-1],]*pi.matrix[,x.vec[i+1]])
out[ndat] <- sample.int(nregimes,1,prob=dfun.no.mp(dat,ndat,out,nregimes,prior.mu,prior.n,prior.alpha,prior.beta)*pi.matrix[out[ndat-1],])
out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
