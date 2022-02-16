
#######################################
# response surface experimental design
#######################################

# data frame for model inputs and outputs
experimental.design <- data.frame(
	focal = c(rep("i",5), rep("j", 5)),
	starting.plants.i = c(3, 2, 1, 1, 1, 2, 1, 0, 0, 0),
	starting.plants.j = c(0, 0, 0, 1, 2, 1, 1, 1, 2, 3),
	time = rep(1, 10),
	stringsAsFactors = FALSE
)

# number of experimental replicates to simulate
nreps <- 10

#####################################
# parameters for the fecundity model
#####################################

# define all parameter values
# paired values are always ordered (i, j)
params <- list(
	gamma = c(0.0, 0.0),
	mu    = c(0.0, 0.0),
	nu    = c(0.0, 0.0),
	r     = c(10.0, 5.0),
	K     = c(1.0, 1.0),
	beta  = c(0.001, 0.001),
	phi   = c(2500, 250),
	alpha_ij = 0.90,
	alpha_ji = 0.70
)

########################################################
# simulate the outcomes the response-surface experiment
########################################################

# set the parameters for the model
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

# add output columns to the experimental design
experimental.outcomes <- experimental.design
newcols <- c(
	"ending.plants",
	"biomass",
	"total.fecundity",
	"per.capita.fecundity"
)
for(cname in newcols){
	for(sp in c("i","j")){
		experimental.outcomes[,paste0(cname,".",sp)] <- NA
	}
}

# run the model for each experimental treatment and save the results
for(rr in 1:nrow(experimental.outcomes)){
	Ni0 <- experimental.outcomes$starting.plants.i[rr]
	Nj0 <- experimental.outcomes$starting.plants.j[rr]
	time <- experimental.outcomes$time[rr]
	
	# starting conditions are seeds in seed bank, plants, and plant biomass
	x0 <- c(
		0,Ni0,Ni0*params$beta[1],
		0,Nj0,Nj0*params$beta[2]
	)
	
	# output order matches the conditions above but after time in first column
	growing.season <- fecundity_dynamics(init=x0, duration=time, step_size=time/1000.)
	
	# name the output columns for convenience
	colnames(growing.season) <- c(
		"Time",
		"seeds.i",
		"plants.i",
		"biomass.i",
		"seeds.j",
		"plants.j",
		"biomass.j"
	)

	# use final plants, final biomass, and conversion rate to estimate total and per capita fecundity of i
	plants_i <- growing.season[nrow(growing.season),"plants.i"]
	biomass_i <- growing.season[nrow(growing.season),"biomass.i"]
	seeds_i <- biomass_i * params$phi[1]
	fecundity_i <- seeds_i / plants_i

	# use final plants, final biomass, and conversion rate to estimate total and per capita fecundity of j
	plants_j <- growing.season[nrow(growing.season),"plants.j"]
	biomass_j <- growing.season[nrow(growing.season),"biomass.j"]
	seeds_j <- biomass_j * params$phi[2]
	fecundity_j <- seeds_j / plants_j

	# fillin the outcomes into the table
	experimental.outcomes[rr,5:12] <- c(
		plants_i,
		plants_j,
		biomass_i,
		biomass_j,
		seeds_i,
		seeds_j,
		fecundity_i,
		fecundity_j
	)
}

##########################################################
# generate random data set based on simulated predictions
##########################################################

# set a random seed for reproducibility
# obtained from random.org (Timestamp: 2022-02-09 20:06:47 UTC)
set.seed(743761)

# use a poisson model to generate count fecundities based on the above simulations
simulated.data <- do.call(rbind,apply(
	experimental.outcomes,
	1,
	function(x,nreps){
		data.frame(
			focal=rep(x["focal"], nreps),
			plants.i=rep(as.integer(x["starting.plants.i"]), nreps),
			plants.j=rep(as.integer(x["starting.plants.j"]), nreps),
			seeds.i=0,
			seeds.j=0,
			time=as.numeric(x["time"]),
			fecundity=rpois(nreps, as.numeric(x[paste0("per.capita.fecundity.",x["focal"])])),
			stringsAsFactors=FALSE
		)
	},
	nreps=nreps
))
rownames(simulated.data) <- 1:nrow(simulated.data)

# save the parameters that were used for simulations
simulated.params <- params
