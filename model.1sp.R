
# we use this to compile the dynamic model
library(odeintr)

# define the ode in C++ format to use odeintr
fecundity.dynamics.sys = '
	// seeds in the seed bank
	dxdt[0] = -gamma * x[0] - mu * x[0];

	// plants that have germinated
	dxdt[1] = gamma * x[0] - nu * x[1];

	// plant biomass
	dxdt[2] = x[2] * r * (1 - x[2] / K) + beta * gamma * x[0] - nu * x[2];
'

odeintr::compile_sys(
	"fecundity_dynamics",
	fecundity.dynamics.sys,
	pars = c("gamma", "mu", "nu", "r", "K", "beta")
)
