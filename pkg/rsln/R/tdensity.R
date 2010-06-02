tdensity <-
function(xx,mu,lambda,alpha) exp(lgamma(0.5*(alpha+1))-lgamma(0.5*alpha))*sqrt(lambda/(alpha*pi))*(1 + lambda*(xx-mu)^2/alpha)^(-(alpha+1)/2)

