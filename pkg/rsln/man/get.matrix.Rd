\name{get.matrix}
\alias{get.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.matrix(mat.name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat.name}{
%%     ~~Describe \code{mat.name} here~~
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
function(mat.name){
if(mat.name=="strong") out <- matrix(c(0.8,0.95,0.2,0.05),2)
if(mat.name=="weak") out <- matrix(c(0.6,0.65,0.4,0.35),2)
if(mat.name=="none") out <- matrix(c(0.5,0.5,0.5,0.5),2)
if(mat.name=="absorb") out <- matrix(c(0.8,0.3,0.2,0.7),2)
out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
