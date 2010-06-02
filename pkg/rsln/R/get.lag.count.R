get.lag.count <-
function(current.x,nregimes) as.matrix(table(factor(current.x[-length(current.x)],levels=1:nregimes),factor(current.x[-1],levels=1:nregimes)))

