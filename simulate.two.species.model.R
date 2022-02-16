
# to use/run the 2sp model
library(here)
source(here('code/R/lib/model.2sp.R'))

#####################
#####################
#####################

# initial conditions for viable seeds in seed bank
N0 <- c(50, 75)

# number of years to simulate
nyears <- 3

# define all parameter values for seed germination phase
# paired values are always ordered (i, j)
params.seed <- list(
	T     = 0.50,
	gamma = c(0.10, 0.15),
	mu    = c(0.01, 0.02),
	nu    = c(0.00, 0.20),
	r     = c(0.00, 0.00),
	K     = c(100.0, 250.0),
	beta  = c(0.02, 0.02),
	phi   = c(5,25),
	alpha_ij = 0.75,
	alpha_ji = 0.50
)

# define all parameter values for plant growth phase
# paired values are always ordered (i, j)
params.vegetative <- list(
	T     = 0.50,
	gamma = c(0.0, 0.01),
	mu    = c(0.00, 0.00),
	nu    = c(0.05, 0.20),
	r     = c(12.00, 12.00),
	K     = c(100.0, 250.0),
	beta  = c(0.02, 0.02),
	phi   = c(5,5),
	alpha_ij = 0.75,
	alpha_ji = 0.50
)

#####################
#####################
#####################

# initial condition including viable seeds, plants, and biomass
x0 <- c(
	N0[1], 0, 0,
	N0[2], 0, 0,
)

for(year in 1:nyears){

	# set the parameters for phase 1
	fecundity_dynamics_set_params(
		gamma_i = params$gamma[1],
		mu_i = params$mu[1],
		nu_i = params$nu[1],
		r_i = params$r[1],
		K_i = params$K[1],
		beta_i = params$beta[1],
		gamma_j = params$gamma[2],
		mu_j = params$mu[2],
		nu_j = params$nu[2],
		r_j = params$r[2],
		K_j = params$K[2],
		beta_j = params$beta[2],
		alpha_ij = params$alpha_ij,
		alpha_ji = params$alpha_ji
	)
	NPBearly <- fecundity_dynamics(x0, params$T[1], params$T[1]/1000.)

	# create germinant biomass

	# set the parameters for phase 2
	fecundity_dynamics_set_params(
		gamma_i = params$gamma[1],
		mu_i = params$mu[1],
		nu_i = params$nu[1],
		r_i = params$r[1],
		K_i = params$K[1],
		beta_i = params$beta[1],
		gamma_j = params$gamma[2],
		mu_j = params$mu[2],
		nu_j = params$nu[2],
		r_j = params$r[2],
		K_j = params$K[2],
		beta_j = params$beta[2],
		alpha_ij = params$alpha_ij,
		alpha_ji = params$alpha_ji
	)
	NPBlate <- fecundity_dynamics(x1, params$T[1], params$T[1]/1000.)


#####################
#####################
#####################

# starting conditions for plant densities
Nis <- seq.int(0,50)
Njs <- seq.int(0,53)

# number of years to simulate
nyears <- 1

#####################
#####################
#####################

# simulate single species density dependent fecundity
fecundity.stats <- expand.grid(Nis, Njs)
colnames(fecundity.stats) <- c("starting.plants.i", "starting.plants.j")
newcols <- c(
	"ending.plants",
	"biomass",
	"total.fecundity",
	"per.capita.fecundity"
)
for(cname in newcols){
	for(sp in c("i","j")){
		fecundity.stats[,paste0(cname,".",sp)] <- NA
	}
}

# to save some time series
bingo <- list()
	
for(rr in 1:nrow(fecundity.stats)){
	Ni0 <- fecundity.stats$starting.plants.i[rr]
	Nj0 <- fecundity.stats$starting.plants.j[rr]
	
	x0 <- c(
		25,Ni0,Ni0*params$beta[1],
		25,Nj0,Nj0*params$beta[2]
	)
	
	growing.season <- fecundity_dynamics(x0, params$T, params$T/1000.)
	
	plants_i <- growing.season[nrow(growing.season),3]
	biomass_i <- growing.season[nrow(growing.season),4]
	seeds_i <- biomass_i * params$phi[1]
	fecundity_i <- seeds_i / plants_i

	plants_j <- growing.season[nrow(growing.season),6]
	biomass_j <- growing.season[nrow(growing.season),7]
	seeds_j <- biomass_j * params$phi[2]
	fecundity_j <- seeds_j / plants_j

	fecundity.stats[rr,3:10] <- c(
		plants_i,
		plants_j,
		biomass_i,
		biomass_j,
		seeds_i,
		seeds_j,
		fecundity_i,
		fecundity_j
	)

	bingo[[paste(Ni0, Nj0, sep=".")]] <- growing.season
}

setEPS()
postscript('./response_surface_2sp.eps', width=6.96, height=7.56)
# par(mfrow=c(2,2), 
par(oma=c(0,0,0,0), mar=c(4.5,5,4,2))

# par(mar = c(4.5, 1.5, 0.5, 3.5), oma = c(0.25, 4.0, 0.0, 0.0))

layout(mat = matrix(
		1:4,
		byrow=TRUE,
		nrow = 2,
		ncol = 2
	),
	heights = c(6,6),
	widths = rep(2,6)
)

par(cex.lab=1.5,cex.axis=1.2)
par(lwd=1.5)


x <- sort(unique(round(fecundity.stats$ending.plants.i,4)))
y <- sort(unique(round(fecundity.stats$ending.plants.j,4)))

# contours of biomass of species i
z <- reshape2::acast(
	fecundity.stats,
	formula = round(ending.plants.i,4) ~ round(ending.plants.j,4),
	value.var="biomass.i"
)
contour(
	x=x,
	y=y,
	z=z,
	labcex=1.25,
	xlim=c(0,50),
	ylim=c(0,50),
	xlab=expression("Plants of species "*italic(i)*" ("*italic(n[i])*")"),
	ylab=expression("Plants of species "*italic(j)*" ("*italic(n[j])*")"),
	# asp=1,
	# method="edge"
	# # nlevels=7
)
title(expression(bold("A")*" Biomass of species "*italic(i)*" ("*italic(B[i])*")"), line=1.75, cex.main=1.5, adj=0.40)

# contours of fecundity of species i
z <- reshape2::acast(
	fecundity.stats,
	formula = round(ending.plants.i,4) ~ round(ending.plants.j,4),
	value.var="per.capita.fecundity.i"
)
contour(
	x=x,
	y=y,
	z=z,
	labcex=1.25,
	xlim=c(0,50),
	ylim=c(0,50),
	xlab=expression("Plants of species "*italic(i)*" ("*italic(n[i])*")"),
	ylab=expression("Plants of species "*italic(j)*" ("*italic(n[j])*")"),
	# method="edge",
	# levels=seq(0,22,2)
)
# title(expression(bold("B")*" Per capita fecundity of species "*italic(i)), line=1.7, cex.main=1.5, adj=0.40)
title(expression(bold("B")), line=1.75, cex.main=1.5, adj=0.00)
title(expression("Per capita fecundity"), line=2.5, cex.main=1.5, adj=0.50)
title(expression("of species "*italic(i)*" ("*italic(F[i])*")"), line=1.00, cex.main=1.5, adj=0.60)

# contours of biomass of species i
z <- reshape2::acast(
	fecundity.stats,
	formula = round(ending.plants.i,4) ~ round(ending.plants.j,4),
	value.var="biomass.j"
)
contour(
	x=x,
	y=y,
	z=z,
	labcex=1.25,
	xlim=c(0,50),
	ylim=c(0,50),
	xlab=expression("Plants of species "*italic(i)*" ("*italic(n[i])*")"),
	ylab=expression("Plants of species "*italic(j)*" ("*italic(n[j])*")"),
	# method="edge"
)
title(expression(bold("C")*" Biomass of species "*italic(j)*" ("*italic(B[j])*")"), line=1.75, cex.main=1.5, adj=0.40)

# contours of fecundity of species i
z <- reshape2::acast(
	fecundity.stats,
	formula = round(ending.plants.i,4) ~ round(ending.plants.j,4),
	value.var="per.capita.fecundity.j"
)
contour(
	x=x,
	y=y,
	z=z,
	labcex=1.25,
	xlim=c(0,50),
	ylim=c(0,50),
	xlab=expression("Plants of species "*italic(i)*" ("*italic(n[i])*")"),
	ylab=expression("Plants of species "*italic(j)*" ("*italic(n[j])*")"),
	# method="simple"
	# levels=c(8,10,12,14,16,20,30)
)
title(expression(bold("D")), line=1.75, cex.main=1.5, adj=0.00)
title(expression("Per capita fecundity"), line=2.5, cex.main=1.5, adj=0.50)
title(expression("of species "*italic(j)*" ("*italic(F[j])*")"), line=1.00, cex.main=1.5, adj=0.60)

# Njs <- c(1, 5, 10, 20)
# for(j in seq_along(Njs)){
# 	hmm <- subset(fecundity.stats, starting.plants.j == Njs[j])

# 	if(j == 1){
# 		plot(
# 			hmm$ending.plants.i,
# 			hmm$biomass.i,
# 			type='l',
# 			col=1,
# 			lty=1,
# 			xlab=expression("Plants of species "*italic(i)*" ("*italic(n[i])*")"),
# 			ylab=expression("Biomass of species "*italic(i)*" ("*italic(B[i])*")"),
# 			xlim=c(0,50),
# 			ylim=c(0,80)
# 		)
# 	}else{
# 		lines(
# 			hmm$ending.plants.i,
# 			hmm$biomass.i,
# 			# type='l',
# 			col=1,
# 			lty=j,
# 			# ylim=c(0,100)
# 		)
# 	}
# }
# title("A", line=-1.7, cex.main=1.5, adj=0.10)

# for(j in seq_along(Njs)){
# 	hmm <- subset(fecundity.stats, starting.plants.j == Njs[j])

# 	if(j == 1){
# 		plot(
# 			hmm$ending.plants.i,
# 			hmm$per.capita.fecundity.i,
# 			type='l',
# 			col=1,
# 			lty=1,
# 			xlab=expression("Plants of species "*italic(i)*" ("*italic(n[i])*")"),
# 			ylab="",
# 			xlim=c(0,50),
# 			ylim=c(0,25)
# 		)
# 	}else{
# 		lines(
# 			hmm$ending.plants.i,
# 			hmm$per.capita.fecundity.i,
# 			# type='l',
# 			col=1,
# 			lty=j,
# 			# ylim=c(0,100)
# 		)
# 	}
# }
# title("B", line=-1.7, cex.main=1.5, adj=0.90)
# mtext(
# 	expression("Per capita fecundity"),
# 	2,
# 	cex=1.2,
# 	line=3.5
# )
# mtext(
# 	expression("of species "*italic(i)*" ("*italic(F[i])*")"),
# 	2,
# 	cex=1.2,
# 	line=2.2
# )

# Nis <- c(1, 5, 10, 20)
# for(i in seq_along(Nis)){
# 	hmm <- subset(fecundity.stats, starting.plants.i == Nis[i])

# 	if(i == 1){
# 		plot(
# 			hmm$ending.plants.j,
# 			hmm$biomass.j,
# 			type='l',
# 			col=1,
# 			lty=1,
# 			xlab=expression("Plants of species "*italic(j)*" ("*italic(n[j])*")"),
# 			ylab=expression("Biomass of species "*italic(j)*" ("*italic(n[j])*")"),
# 			xlim=c(0,50),
# 			ylim=c(0,80)
# 		)
# 	}else{
# 		lines(
# 			hmm$ending.plants.j,
# 			hmm$biomass.j,			
# 			# type='l',
# 			col=1,
# 			lty=i,
# 			# ylim=c(0,100)
# 		)
# 	}
# }
# title("C", line=-1.7, cex.main=1.5, adj=0.10)

# for(i in seq_along(Nis)){
# 	hmm <- subset(fecundity.stats, starting.plants.i == Nis[i])

# 	if(i == 1){
# 		plot(
# 			hmm$ending.plants.j,
# 			hmm$per.capita.fecundity.j,			
# 			type='l',
# 			col=1,
# 			lty=1,
# 			xlab=expression("Plants of species "*italic(j)*" ("*italic(n[j])*")"),
# 			ylab="",
# 			xlim=c(0,50),
# 			ylim=c(0,150)
# 		)
# 	}else{
# 		lines(
# 			hmm$ending.plants.j,
# 			hmm$per.capita.fecundity.j,			
# 			# type='l',
# 			col=1,
# 			lty=i,
# 			# ylim=c(0,100)
# 		)
# 	}
# }
# title("D", line=-1.7, cex.main=1.5, adj=0.90)
# mtext(
# 	expression("Per capita fecundity"),
# 	2,
# 	cex=1.2,
# 	line=3.5
# )
# mtext(
# 	expression("of species "*italic(j)*" ("*italic(F[j])*")"),
# 	2,
# 	cex=1.2,
# 	line=2.2
# )

dev.off()
