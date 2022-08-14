# About
This repository contains the Master's thesis and the R code pertaing to it. The thesis was written as part of the Master's in Statistical Science and Data Science at Leiden University between February and August $2022$.

This research compares three commonly used matching methods; Nearest neighbour matching with replacement and without replacement, and propensity score matching using a caliper (PSM). Monte Carlo simulations were performed to assess the strengths and weaknesses of the methods in estimating the true treatment effect under different conditions.


## Thesis: Match Making
Pdf file of master's thesis


## [R_functions](/matching_simulation_functions.R)

Functions for performing the simulations

## [Uniform_simulations_1](https://github.com/laura-ruth/LJS_thesis/blob/main/simulations1.Rmd)

MC simulations for uniformly distributed data comparing 3 matching methods. Parameters are sample size, ratio, caliper size, and treatment effect size. 

## [Simulations_ratio_1_to_1](https://github.com/laura-ruth/LJS_thesis/blob/main/unif_large_n_ratio1.Rmd)

MC simulations comparing ther 3 matching methods when the ratio is 1:1 for sample sizes from 10 to 1010. 

## [Uniform_simulations_2_overlap](/simulations2_overlap_uniform.Rmd)

MC simulations on overlap for the uniform distribution

## [Interaction_simulations](https://github.com/laura-ruth/LJS_thesis/blob/main/Interaction_Simulations.Rmd)

MC simulations looking at the effect of different levels of interaction on the estimates produced by the 3 matching methods.

## [Propensity_simulations](https://github.com/laura-ruth/LJS_thesis/blob/main/propensity%20_simulations.Rmd)

Investigations into the propensity score as implemented in R.
