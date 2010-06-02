rdirichlet <-
function(alpha.vec) {
gamma.vec <- sapply(alpha.vec, function(xx) rgamma(1,xx))
gamma.vec/sum(gamma.vec)
}

