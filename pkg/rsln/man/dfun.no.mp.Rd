\name{dfun.no.mp}
\alias{dfun.no.mp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dfun.no.mp(dat, index, x.vec, nregimes, prior.mu, prior.n, prior.alpha, prior.beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
}
  \item{x.vec}{
%%     ~~Describe \code{x.vec} here~~
}
  \item{nregimes}{
%%     ~~Describe \code{nregimes} here~~
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
