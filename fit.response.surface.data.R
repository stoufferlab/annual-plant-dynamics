
# libraries required to run the code
library(bbmle)
library(dotwhisker)
library(dplyr)
library(odeintr)

# to use/run the 2sp model
source('lib/model.2sp.R')

# create the simulated dataset
# see file below for response surface "exerimental design"
source('lib/generate.response.surface.data.R')

# load the functions used for model fitting
source('lib/model.fitting.functions.R')

#################################
# example model fitting approach
#################################

# first estimate the best-fit model parameters for i only data

# starting guesses for i only parameters
params <- list(
	log_r_i = log(10),
	log_phi_i = log(max(subset(simulated.data, focal=="i")$fecundity))
)

# use mle2 to infer unknown values for i only data
# NOTE: unidentifiable parameters are held at fixed values
model.fit.i <- bbmle::mle2(
	minuslogl = fecundity.model.NLL,
	start = unlist(params),
	fixed = list(
	 	gamma_i = 0,
		gamma_j = 0,
		mu_i = 0,
		mu_j = 0,
		nu_i = 0,
		nu_j = 0,
		beta_i = 0.001,
		beta_j = 0.001,
		K_i = 1,
		K_j = 1,
		log_r_j = 0,
		log_phi_j = 0,
		alpha_ij = 0,
		alpha_ji = 0
	),
	data=subset(simulated.data, plants.j == 0),
	vecpar = TRUE,
	control = list(
		maxit=100000
	),
	skip.hessian=FALSE
)

# second estimate the best-fit model parameters for j only data

# starting guesses for j only parameters
params$log_r_j <- log(10)
params$log_phi_j = log(max(subset(simulated.data, focal=="j")$fecundity))

# use mle2 to infer unknown values for j only data
# NOTE: unidentifiable parameters are held at fixed values
model.fit.j <- bbmle::mle2(
	minuslogl = fecundity.model.NLL,
	start = unlist(params),
	fixed = list(
	 	gamma_i = 0,
		gamma_j = 0,
		mu_i = 0,
		mu_j = 0,
		nu_i = 0,
		nu_j = 0,
		beta_i = 0.001,
		beta_j = 0.001,
		K_i = 1,
		K_j = 1,
		log_r_i = 0,
		log_phi_i = 0,
		alpha_ij = 0,
		alpha_ji = 0
	),
	data=subset(simulated.data, plants.i == 0),
	vecpar = TRUE,
	control = list(
		maxit=100000
	),
	skip.hessian=FALSE
)

# lastly we will use the estimated parameters above as starting values when fitting ALL data
params$log_r_i <- as.numeric(coef(model.fit.i)["log_r_i"])
params$log_r_j <- as.numeric(coef(model.fit.j)["log_r_j"])
params$log_phi_i <- as.numeric(coef(model.fit.i)["log_phi_i"])
params$log_phi_j <- as.numeric(coef(model.fit.j)["log_phi_j"])

# assume no heterospecific interactions as starting guess
# (admittedly, this is likely a poor guess)
params$alpha_ij <- 0
params$alpha_ji <- 0

# use mle2 to estimate ALL best-fit model parameters and associated SEs
# NOTE: unidentifiable parameters are held at fixed values
model.fit <- bbmle::mle2(
	minuslogl = fecundity.model.NLL,
	start = unlist(params),
	fixed = list(
	 	gamma_i = 0,
		gamma_j = 0,
		mu_i = 0,
		mu_j = 0,
		nu_i = 0,
		nu_j = 0,
		K_i = 1,
		K_j = 1,
		beta_i = 0.001,
		beta_j = 0.001
	),
	data=simulated.data,
	vecpar = TRUE,
	control = list(
		maxit=100000
	),
	skip.hessian=FALSE
)

# print out the coefficient summary table
print(summary(model.fit))

# evaluate a model with the parameters used for the simulations
known.params <- list()
known.params$log_r_i <- log(simulated.params$r[1])
known.params$log_r_j <- log(simulated.params$r[2])
known.params$log_phi_i <- log(simulated.params$phi[1])
known.params$log_phi_j <- log(simulated.params$phi[2])
known.params$alpha_ij <- simulated.params$alpha_ij
known.params$alpha_ji <- simulated.params$alpha_ji

# use mle2 to calculate loglikelihood, but mostly to generate a coefficient table
model.known <- bbmle::mle2(
	minuslogl = fecundity.model.NLL,
	start = unlist(known.params),
	fixed = list(
	 	gamma_i = 0,
		gamma_j = 0,
		mu_i = 0,
		mu_j = 0,
		nu_i = 0,
		nu_j = 0,
		K_i = 1,
		K_j = 1,
		beta_i = 0.001,
		beta_j = 0.001
	),
	data=simulated.data,
	vecpar = TRUE,
	eval.only=TRUE
)

# plot known parameters versus inferred parameters
dwplot(list(model.known, model.fit)) %>%
relabel_predictors(
	c(
		log_r_i = "ln(r_i)",
		log_r_j = "ln(r_j)",
		log_phi_i = "ln(phi_i)",
		log_phi_j = "ln(phi_j)"
	)
) +
scale_colour_discrete(
    name = "Parameters",
    labels = c("Inferred", "Known")
)
