# Annual-Plant Population Dynamics

This repository contains code to support the unpublished study *A critical examination of models of annual-plant population dynamics and density-dependent fecundity* by Daniel B. Stouffer (daniel.stouffer@canterbury.ac.nz).

## Organization
The repository provides two example scripts (`simulate.two.species.model.R` and `fit.response.surface.data.R`) and additional scripts in [lib](lib/) which provide utility functions used therein.

Succinctly, the code in `simulate.two.species.model.R` allows the user to simulate continuous-time population dynamics of two co-occurring annual-plant species using the model developed in the publication. The code in `fit.response.surface.data.R` provides an example for how that same model could be fit to simulated data designed to mirror a common two-species response--surface experimental design.

## Warranty
All code is provided "as is" and without warranty. If you know of more efficient or elegant ways of doing anything (of which there are likely many), please let me know.
