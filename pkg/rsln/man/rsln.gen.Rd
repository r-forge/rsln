\name{rsln.gen}
\alias{rsln.gen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rsln.gen(ndata, nregimes, mean.vec = rep(0, nregimes), prec.vec = rep(1, nregimes), trans.mat = matrix(1/nregimes, nregimes, nregimes))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ndata}{
%%     ~~Describe \code{ndata} here~~
}
  \item{nregimes}{
%%     ~~Describe \code{nregimes} here~~
}
  \item{mean.vec}{
%%     ~~Describe \code{mean.vec} here~~
}
  \item{prec.vec}{
%%     ~~Describe \code{prec.vec} here~~
}
  \item{trans.mat}{
%%     ~~Describe \code{trans.mat} here~~
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
function(ndata,nregimes,mean.vec=rep(0,nregimes),prec.vec=rep(1,nregimes),trans.mat=matrix(1/nregimes,nregimes,nregimes)){
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
