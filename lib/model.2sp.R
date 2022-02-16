
##################################
# model simulation functions
##################################

# we use this to compile the dynamic model
library(odeintr)

# define the ode in C++ format to use odeintr
fecundity.dynamics.sys = '
	// seeds in the seed bank species i
	dxdt[0] = -gamma_i * x[0] - mu_i * x[0];

	// plants that have germinated species i
	dxdt[1] = gamma_i * x[0] - nu_i * x[1];

	// plant biomass species i
	dxdt[2] = x[2] * r_i * (1 - (x[2] + alpha_ij * x[5]) / K_i)
	        + beta_i * gamma_i * x[0]
	        - nu_i * x[2];

	// seeds in the seed bank species j
	dxdt[3] = -gamma_j * x[3] - mu_j * x[3];

	// plants that have germinated species j
	dxdt[4] = gamma_j * x[3] - nu_j * x[4];

	// plant biomass species j
	dxdt[5] = x[5] * r_j * (1 - (alpha_ji * x[2] + x[5]) / K_j)
	        + beta_j * gamma_j * x[3] - nu_j * x[5];
'

# compile the model into the local environment
odeintr::compile_sys(
	"fecundity_dynamics",
	fecundity.dynamics.sys,
	pars = c(
		"gamma_i",
		"gamma_j",
		"mu_i",
		"mu_j",
		"nu_i",
		"nu_j",
		"r_i",
		"r_j",
		"K_i",
		"K_j",
		"beta_i",
		"beta_j",
		"alpha_ij",
		"alpha_ji"
	),
	const = TRUE,
	method = 'bsd'
)
