
##################################
# model fitting functions
##################################

# produce two-species model prediction given:
# 1. model parameters
# 2. a set of initial conditions for all state variables
fecundity.model.predict = function(params, plants.i, plants.j, seeds.i, seeds.j, time, focal, verbose=FALSE){

	# set the parameters in the ode solver
	fecundity_dynamics_set_params(
		gamma_i = params["gamma_i"],
		mu_i = params["mu_i"],
		nu_i = params["nu_i"],
		r_i = exp(params["log_r_i"]),
		K_i = params["K_i"],
		beta_i = params["beta_i"],
		gamma_j = params["gamma_j"],
		mu_j = params["mu_j"],
		nu_j = params["nu_j"],
		r_j = exp(params["log_r_j"]),
		K_j = params["K_j"],
		beta_j = params["beta_j"],
		alpha_ij = params["alpha_ij"],
		alpha_ji = params["alpha_ji"]
	)

	# for tracing purposes
	if(verbose>1){
		print(params)
		message("predicting")
		flush(stdout())
	}

	# container for the predicted outputs
	predicted.fecundity <- numeric(length(plants.i))

	# iterate over all observations
	for(i in seq.int(length(plants.i))){
		# for tracing purposes
		if(verbose>2){
			message(i)
			flush(stdout())
		}

		# starting conditions are viable seeds in seed bank, plants, and plant biomass
		x0 <- c(
			seeds.i[i], plants.i[i], plants.i[i] * params["beta_i"],
			seeds.j[i], plants.j[i], plants.j[i] * params["beta_j"]
		)
	
		# output order matches the conditions above but with time elapsed in first column
		growing.season <- fecundity_dynamics(x0, time[i], time[i]/1000.)
		colnames(growing.season) <- c(
			"time.elapsed",
			"seeds.i",
			"plants.i",
			"biomass.i",
			"seeds.j",
			"plants.j",
			"biomass.j"
		)

		# we only want the final values to make our prediction
		growing.season <- growing.season[nrow(growing.season),]

		# convert biomass to per capita fecundity
		growing.season$fecundity.i <- exp(params["log_phi_i"]) * growing.season$biomass.i / growing.season$plants.i
		growing.season$fecundity.j <- exp(params["log_phi_j"]) * growing.season$biomass.j / growing.season$plants.j

		# select the fecundity corresponding to the focal species
		predicted.fecundity[i] <- growing.season[,paste0("fecundity.",focal[i])]
	}

	# for tracing purposes
	if(verbose>1){
		message("predicted")
		flush(stdout())
	}

	return(predicted.fecundity)
}

# calculate the negative loglikelihood of a set of observations given:
# 1. model parameters
# 2. a set of initial conditions for all state variables
# 3. observed fecundities
fecundity.model.NLL = function(params, plants.i, plants.j, seeds.i, seeds.j, time, focal, fecundity, verbose=0){
	# calculate the vector of predicted values using function above
	predicted.fecundity <- fecundity.model.predict(
		params,
		plants.i = plants.i,
		plants.j = plants.j,
		seeds.i = seeds.i,
		seeds.j = seeds.j,
		time = time,
		focal = focal,
		verbose = verbose
	)
	
	# anywhere in parameter space that is non-biolgical or uninformative should be avoided
	# otherwise we treat observed fecundities as Poisson observations to calculate the log-likelihood
	if(!all(predicted.fecundity > 0) || any(!is.finite(predicted.fecundity))){
		return(Inf)
	}else{
		nll <- -sum(dpois(fecundity, predicted.fecundity, log=TRUE))
	}

	# for tracing purposes
	if(verbose>0){
		print(nll)
		flush(stdout())
	}

	return(nll)
}

# define the order of parameters ; this is a requirement of mle2 to use a parameter vector (like optim)
parnames(fecundity.model.NLL) <- c(
	"gamma_i",
	"gamma_j",
	"mu_i",
	"mu_j",
	"nu_i",
	"nu_j",
	"log_r_i",
	"log_r_j",
	"K_i",
	"K_j",
	"beta_i",
	"beta_j",
	"log_phi_i",
	"log_phi_j",
	"alpha_ij",
	"alpha_ji"
)
