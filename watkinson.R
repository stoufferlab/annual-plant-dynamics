
# library(RColorBrewer)

# libraries used below
library(odeintr)

# define the ode in C++ format to use odeintr
fecundity.dynamics.sys = '
	// seeds in the seed bank
	dxdt[0] = -g * x[0] - m_s * x[0];

	// plants that have germinated
	dxdt[1] = g * x[0] - m_p * x[1];

	// plant biomass
	dxdt[2] = x[2] * r * (1 - a * x[2] / K) - m_p * x[2] + b * g * x[0];
'

odeintr::compile_sys(
	"fecundity_dynamics",
	fecundity.dynamics.sys,
	pars = c("g", "m_s", "m_p", "r", "K", "a", "b")
	# method = "bsd"
)

#####################
#####################
#####################

# define all parameter values
params_early <- list(
	g    = 1.5,
	m_s  = 0.75,
	m_p  = 0.0,
	r    = 0,
	K    = 100,
	a    = 1,
	b    = .02
)

params_middle <- list(
	g    = 0.0,
	m_s  = 0.0,
	m_p  = 0.0,
	r    = 12,
	K    = 100,
	a    = 1,
	b    = .02
)

# params_late <- list(
# 	g    = 0.0,
# 	m_s  = 0.0,
# 	m_p  = 0.0,
# 	r    = 0.0,
# 	K    = 100,
# 	a    = 1,
# 	b    = .02
# )

conversion <- 0.85

# length of different time periods
year.length <- 1
T1 <- year.length/2
T2 <- year.length/2
# T3 <- year.length/4

# starting conditions
x0 <- c(100,0,0)

# number of years to simulate
nyears <- 2

SPtotal <- data.frame()
seed.stats <- matrix(NA,nrow=nyears,ncol=4)
fecundity.stats <- matrix(NA,nrow=nyears,ncol=4)
for(year in seq.int(nyears)){
	# simulate seed / plant germination and seed mortality for the first time period
	with(params_early,
		fecundity_dynamics_set_params(
			g = g,
			m_s	= m_s,
			m_p = m_p,
			r = r,
			K = K,
			a = a,
			b = b
		)
	)

	SPearly <- fecundity_dynamics(x0, T1, T1/1000.)
	SPearly$Time <- SPearly$Time + year.length*(year-1)

	# simulate seed / plant germination, seed / plant mortality, and plant growth for the second time period
	x1 <- unlist(tail(SPearly[,2:ncol(SPearly)], 1))

	# use x0 and x1 to determine the seeds that germinated vs those that died
	seed.stats[year,] <- c(
		T1 + year.length*(year-1),
		x0[1],
		x1[2],
		x0[1]-(x1[1]+x1[2])
	)

	with(params_middle,
		fecundity_dynamics_set_params(
			g = g,
			m_s	= m_s,
			m_p = m_p,
			r = r,
			K = K,
			a = a,
			b = b
		)
	)

	SPmid <- fecundity_dynamics(x1, T2, T2/100.)
	SPmid$Time <- SPmid$Time + T1 + year.length*(year-1)

	# growing season ends: proportion of biomass is converted to seeds and all plants undergo senescence 
	x2 <- unlist(tail(SPmid[,2:ncol(SPmid)], 1))

	# rapid transformation step
	x2prime <- x2
	x2prime[1] <- x2[1] + conversion * x2[3]
	x2prime[2] <- 0
	x2prime[3] <- 0

	fecundity.stats[year,] <- c(
		T1 + T2 + year.length*(year-1),
		x1[1],
		x1[1]-x2[1],
		conversion*x2[3]
	)

	# with(params_late,
	# 	fecundity_dynamics_set_params(
	# 		g = g,
	# 		m_s	= m_s,
	# 		m_p = m_p,
	# 		r = r,
	# 		K = K,
	# 		a = a,
	# 		b = b
	# 	)
	# )

	# SPlate <- fecundity_dynamics(x2prime, T3, T3/100.)
	# SPlate$Time <- SPlate$Time + T1 + T2 + year.length*(year-1)

	# stitch together and plot

	# SPtotal <- rbind(SPtotal, rbind(SPearly, SPmid, SPlate))
	SPtotal <- rbind(SPtotal, rbind(SPearly, SPmid))

	x0 <- x2prime
}

buffer <- 0.00
# graphics.off()
# dev.new()
setEPS()
postscript('./watkinson.eps', width=9, height=6.5)
# par(mfrow=c(2,2), 
par(oma=c(0,0,0,0), mar=c(4,5,1,1))

layout(mat = matrix(
		1:4,
		byrow=TRUE,
		nrow = 2,
		ncol = 2
	),
	heights = rep(2,6),
	widths = rep(2,6)
)

par(cex.lab=1.45,cex.axis=1.2)
par(lwd=1.5)

plot(
	SPtotal[,1],
	SPtotal[,2],
	type='n',
	col=1,
	lty=3,
	xlab="",
	ylab="Viable seeds in the seed bank",
	ylim=c(0,150)
)

# indicate germination vs mortality
library(shape)
for(year in seq.int(nyears)){
	# T1
	Arrows(
		ifelse(year==1,seed.stats[year,1]-T1,seed.stats[year,1]-T1+buffer),
		seed.stats[year,2],
		seed.stats[year,1],
		seed.stats[year,2],
		code=2,
		arr.type="triangle",
		arr.adj=1,
		arr.length=0.2,
		arr.width=0.175,
		arr.col=1,
		lcol=1,
		lty=1
	)
	# germination
	Arrows(
		seed.stats[year,1],
		seed.stats[year,2],
		seed.stats[year,1],
		seed.stats[year,2]-seed.stats[year,3],
		# segment=FALSE,
		code=2,
		arr.type="triangle",
		arr.adj=1,
		arr.length=0.2,
		arr.width=0.175,
		arr.col=1,
		lcol=1,
		lty=1
	)
	# mortality
	Arrows(
		seed.stats[year,1],
		seed.stats[year,2]-seed.stats[year,3],
		seed.stats[year,1],
		seed.stats[year,2]-seed.stats[year,3]-seed.stats[year,4],
		# segment=FALSE,
		code=2,
		arr.type="triangle",
		arr.adj=1,
		arr.length=0.2,
		arr.width=0.175,
		arr.col=1,
		lcol=1,
		lty=1
	)
	# T2
	Arrows(
		fecundity.stats[year,1]-T2,
		seed.stats[year,2]-seed.stats[year,3]-seed.stats[year,4],
		fecundity.stats[year,1]-buffer,
		fecundity.stats[year,2],
		code=2,
		arr.type="triangle",
		arr.adj=1,
		arr.length=0.2,
		arr.width=0.175,
		arr.col=1,
		lcol=1,
		lty=1
	)

	# # germination
	# Arrows(
	# 	fecundity.stats[year,1]-buffer,
	# 	fecundity.stats[year,2],
	# 	fecundity.stats[year,1]-buffer,
	# 	fecundity.stats[year,2]-fecundity.stats[year,3],
	# 	# segment=FALSE,
	# 	code=2,
	# 	arr.type="triangle",
	# 	arr.adj=1,
	# 	arr.length=0.2,
	# 	arr.width=0.175,
	# 	arr.col=1,
	# 	lcol=1,
	# 	lty=1
	# )

	# fecundity
	Arrows(
		fecundity.stats[year,1]+.01,
		fecundity.stats[year,2]-fecundity.stats[year,3],
		fecundity.stats[year,1]+.01,
		fecundity.stats[year,2]-fecundity.stats[year,3]+fecundity.stats[year,4],
		# segment=FALSE,
		code=2,
		arr.type="triangle",
		arr.adj=1,
		arr.length=0.2,
		arr.width=0.175,
		arr.col=1,
		lcol=1,
		lty=1
	)


}
Arrows(
	fecundity.stats[year,1],
	fecundity.stats[year,2]-fecundity.stats[year,3]+fecundity.stats[year,4],
	fecundity.stats[year,1]+T3,
	fecundity.stats[year,2]-fecundity.stats[year,3]+fecundity.stats[year,4],
		code=2,
		arr.type="triangle",
		arr.adj=1,
		arr.length=0.2,
		arr.width=0.175,
		arr.col=1,
		lcol=1,
		lty=1
	)
# lines(SPtotal[,1], SPtotal[,2], col=1, lty=2)

# seeds (continuous)
plot(
	SPtotal[,1],
	SPtotal[,2],
	type='l',
	col=1,
	lty=2,
	xlab="",
	ylab="Viable seeds in the seed bank",
	ylim=c(0,150)
)

# plants
plot(
	SPtotal[,1],
	SPtotal[,3],
	type='l',
	col=1,
	lty=2,
	xlab="Years",
	ylab="Seed-producing plants",
	ylim=c(0,50)
)

# plant biomass
plot(
	x=SPtotal[,1],
	y=SPtotal[,4],
	# y=ifelse(SPtotal[,3]!=0,SPtotal[,4]/SPtotal[,3],0),
	type='l',
	col=1,
	lty=3,
	xlab="Years",
	ylab="Total plant biomass",
	log='',
	ylim=c(0,100)
)

dev.off()

# plot(SPtotal[,1], ifelse(SPtotal[,3]!=0,SPtotal[,4]/SPtotal[,3],0), type='l', col=1, lty=4, xlab="Years", ylab="Biomass per plant")

# setEPS()
# postscript('./fecundity_twospecies.eps', width=6.96, height=6.6)

# layout(mat = matrix(
# 		1:4,
# 		byrow=TRUE,
# 		nrow = 2,
# 		ncol = 2
# 	),
# 	heights = rep(2,6),
# 	widths = rep(2,6)
# )

# par(mar = c(4.5, 1.5, 0.5, 3.5), oma = c(0.25, 4.0, 0.0, 0.0))

# # set parameters within ode solver
# f_i <- 0
# f_j <- 0
# a_ii <- 0
# a_ij <- 0
# a_ji <- 0
# a_jj <- 0

# b_ii <- 1/15
# b_jj <- 1/15

# p_i <- -(log(1-5*b_ii))/b_ii
# p_j <- -(log(1-5*b_jj))/b_jj

# # tweak the interspecific effects

# b_ji <- 0.75/15
# b_ij <- 0.5/15

# n_i <- seq(1, 100, length.out=1000)
# counter <- 1
# for(n_j in c(0, 5, 10, 20)){
# 	Sfinal <- data.frame() #matrix(NA, nrow=length(n_i), ncol=2)
# 	for(i in seq_along(n_i)){
# 		prod_2sp_set_params(
# 			ni=n_i[i],
# 			nj=n_j,
# 			pi=p_i,
# 			pj=p_j,
# 			fi=f_i,
# 			fj=f_j,
# 			b_ii=b_ii,
# 			b_ij=b_ij,
# 			b_ji=b_ji,
# 			b_jj=b_jj,
# 			a_ii=a_ii,
# 			a_ij=a_ij,
# 			a_ji=a_ji,
# 			a_jj=a_jj
# 		)
# 		S0 <- c(0,0)
# 		St <- prod_2sp(S0, 1, 1/1000.)
# 		Sfinal <- rbind(Sfinal, St[nrow(St),2:3])
# 	}
# 	if(counter == 1){
# 		plot(
# 			Sfinal[,1]/n_i ~ n_i,
# 			type='l',
# 			log='x',
# 			ylim=c(0,8),
# 			xlim=c(1,20),
# 			lty=1,
# 			lwd=1.5,
# 			xaxs='i',
# 			yaxs='i',
# 			las=1,
# 			# axes=FALSE,
# 			bty='L',
# 			xlab='',
# 			ylab=''
# 			# asp=1.0
# 		)
# 	}else{
# 		lines(n_i, Sfinal[,1]/n_i, lty=counter, lwd=1.5)
# 	}
# 	counter <- counter + 1
# }

# title("A", line=-1.7, cex.main=1.5, adj=0.10)

# mtext(
# 	expression(italic(F[i*plain(',')*t])),
# 	2,
# 	line=2.9,
# 	cex=1.3,
# 	las=1
# 	# adj=-0.35
# )

# n_j <- seq(1, 100, length.out=1000)
# counter <- 1
# for(n_i in c(0, 5, 10, 20)){
# 	Sfinal <- data.frame() #matrix(NA, nrow=length(n_i), ncol=2)
# 	for(i in seq_along(n_j)){
# 		prod_2sp_set_params(
# 			ni=n_i,
# 			nj=n_j[i],
# 			pi=p_i,
# 			pj=p_j,
# 			fi=f_i,
# 			fj=f_j,
# 			b_ii=b_ii,
# 			b_ij=b_ij,
# 			b_ji=b_ji,
# 			b_jj=b_jj,
# 			a_ii=a_ii,
# 			a_ij=a_ij,
# 			a_ji=a_ji,
# 			a_jj=a_jj
# 		)
# 		S0 <- c(0,0)
# 		St <- prod_2sp(S0, 1, 1/1000.)
# 		Sfinal <- rbind(Sfinal, St[nrow(St),2:3])
# 	}
# 	if(counter == 1){
# 		plot(
# 			Sfinal[,2]/n_j ~ n_j,
# 			type='l',
# 			log='x',
# 			ylim=c(0,8),
# 			xlim=c(1,20),
# 			lty=1,
# 			lwd=1.5,
# 			xaxs='i',
# 			yaxs='i',
# 			las=1,
# 			# axes=FALSE,
# 			bty='L',
# 			xlab='',
# 			ylab=''
# 			# asp=1.0
# 		)
# 	}else{
# 		lines(n_j, Sfinal[,2]/n_j, lty=counter, lwd=1.5)
# 	}
# 	counter <- counter + 1
# }

# title("B", line=-1.7, cex.main=1.5, adj=0.10)

# mtext(
# 	expression(italic(F[j*plain(',')*t])),
# 	2,
# 	line=2.9,
# 	cex=1.3,
# 	las=1
# 	# adj=-0.35
# )

# # tweak the interspecific effects

# f_j <- 2.5/15

# n_i <- seq(1, 100, length.out=1000)
# counter <- 1
# for(n_j in c(0, 5, 10, 20)){
# 	Sfinal <- data.frame() #matrix(NA, nrow=length(n_i), ncol=2)
# 	for(i in seq_along(n_i)){
# 		prod_2sp_set_params(
# 			ni=n_i[i],
# 			nj=n_j,
# 			pi=p_i,
# 			pj=p_j,
# 			fi=f_i,
# 			fj=f_j,
# 			b_ii=b_ii,
# 			b_ij=b_ij,
# 			b_ji=b_ji,
# 			b_jj=b_jj,
# 			a_ii=a_ii,
# 			a_ij=a_ij,
# 			a_ji=a_ji,
# 			a_jj=a_jj
# 		)
# 		S0 <- c(0,0)
# 		St <- prod_2sp(S0, 1, 1/1000.)
# 		Sfinal <- rbind(Sfinal, St[nrow(St),2:3])
# 	}
# 	if(counter == 1){
# 		plot(
# 			Sfinal[,1]/n_i ~ n_i,
# 			type='l',
# 			log='x',
# 			ylim=c(0,8),
# 			xlim=c(1,20),
# 			lty=1,
# 			lwd=1.5,
# 			xaxs='i',
# 			yaxs='i',
# 			las=1,
# 			# axes=FALSE,
# 			bty='L',
# 			xlab='',
# 			ylab=''
# 			# asp=1.0
# 		)
# 	}else{
# 		lines(n_i, Sfinal[,1]/n_i, lty=counter, lwd=1.5)
# 	}
# 	counter <- counter + 1
# }

# title("C", line=-1.7, cex.main=1.5, adj=0.10)

# mtext(
# 	expression(italic(n[i*plain(',')*t])),
# 	1,
# 	line=2.9,
# 	cex=1.3,
# 	# adj=-0.35
# )

# mtext(
# 	expression(italic(F[i*plain(',')*t])),
# 	2,
# 	line=2.9,
# 	cex=1.3,
# 	las=1
# 	# adj=-0.35
# )

# n_j <- seq(1, 100, length.out=1000)
# counter <- 1
# for(n_i in c(0, 5, 10, 20)){
# 	Sfinal <- data.frame() #matrix(NA, nrow=length(n_i), ncol=2)
# 	for(i in seq_along(n_j)){
# 		prod_2sp_set_params(
# 			ni=n_i,
# 			nj=n_j[i],
# 			pi=p_i,
# 			pj=p_j,
# 			fi=f_i,
# 			fj=f_j,
# 			b_ii=b_ii,
# 			b_ij=b_ij,
# 			b_ji=b_ji,
# 			b_jj=b_jj,
# 			a_ii=a_ii,
# 			a_ij=a_ij,
# 			a_ji=a_ji,
# 			a_jj=a_jj
# 		)
# 		S0 <- c(0,0)
# 		St <- prod_2sp(S0, 1, 1/1000.)
# 		Sfinal <- rbind(Sfinal, St[nrow(St),2:3])
# 	}
# 	if(counter == 1){
# 		plot(
# 			Sfinal[,2]/n_j ~ n_j,
# 			type='l',
# 			log='x',
# 			ylim=c(0,8),
# 			xlim=c(1,20),
# 			lty=1,
# 			lwd=1.5,
# 			xaxs='i',
# 			yaxs='i',
# 			las=1,
# 			# axes=FALSE,
# 			bty='L',
# 			xlab='',
# 			ylab=''
# 			# asp=1.0
# 		)
# 	}else{
# 		lines(n_j, Sfinal[,2]/n_j, lty=counter, lwd=1.5)
# 	}
# 	counter <- counter + 1
# }

# title("D", line=-1.7, cex.main=1.5, adj=0.10)

# mtext(
# 	expression(italic(n[j*plain(',')*t])),
# 	1,
# 	line=2.9,
# 	cex=1.3,
# 	# adj=-0.35
# )

# mtext(
# 	expression(italic(F[j*plain(',')*t])),
# 	2,
# 	line=2.9,
# 	cex=1.3,
# 	las=1
# 	# adj=-0.35
# )

# dev.off()
