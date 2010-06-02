get.stationary.dist <-
function(mat) {
e.mat <- eigen(t(mat))
e.vec <- e.mat$vectors[,zapsmall(e.mat$values)==1]
e.vec/sum(e.vec)
}

