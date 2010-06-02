get.matrix <-
function(mat.name){
if(mat.name=="strong") out <- matrix(c(0.8,0.95,0.2,0.05),2)
if(mat.name=="weak") out <- matrix(c(0.6,0.65,0.4,0.35),2)
if(mat.name=="none") out <- matrix(c(0.5,0.5,0.5,0.5),2)
if(mat.name=="absorb") out <- matrix(c(0.8,0.3,0.2,0.7),2)
out
}

