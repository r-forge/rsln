init.values <- function(){
	nregimes <<- 2
	data.values <<- rsln.gen(100,nregimes)$data.values
	burnin<<-1000
	ndraws<<-1000
	prior.mu<<-0
	prior.n<<-0.001
	prior.alpha<<-1000
	prior.beta<<-0.001
	prior.dirichlet<<-matrix(1,nregimes,nregimes)
	pi.start<<-matrix(1/nregimes,nregimes,nregimes)
}

#init.values()
data.rsln <- rsln.gen(100,3,c(1,-1,3))
system.time(my.fit <- rsln.fit(data.rsln$data.values,2))