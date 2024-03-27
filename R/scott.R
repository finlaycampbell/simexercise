## DATA PREP

## load packages
pacman::p_load(XML, here, data.table, countrycode, readxl, glue, tidyverse, rio,
               deSolve, magrittr)

## source functions
source(here("R/functions.R"))

## get various scenarios defining intervention effects
scenarios <- get_scenarios()

## load country data
cdat <- import(here("data/cdat.rda"))



## PLOT OF A SINGLE SCENARIO

## define a se of parameters
pars <- get_pars(
  R0 = 2.5,
  generation_time = 7,
  incubation_period = 3,
  inf_E_over_I = 0.25,
  ## ifr calculate from meta-analysis
  ifr = age_to_ifr(c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5,
                     42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 80)),
  ## given you go to hospital, what proportion dies?
  hosp_mortality = 1/c(rep(20, 4), rep(15, 4), rep(10, 4), rep(5, 4)),
  ## calculate proportion hospitalised from IFR and proportion in hospital that die
  prop_hosp = ifr/hosp_mortality,
  ## proportion of deaths that going to hospital averts
  hosp_protection = 0.75,
  hosp_duration = c(rep(7, 4), rep(14, 6), rep(21, 6)),
  hosp_capacity = 0.0025,
  comm_mortality = 0,
  vax_rate = 0.002,
  vax_infectious = 0.3,
  vax_infection = 0.5,
  vax_hosp = 0.5,
  vax_death = 0.8,
  vax_prioritised = TRUE,
  hosp_prioritised = TRUE,
  frac_symp = 0.8,
  max_day = 365,
  scenario_type = "Baseline",
  iso3 = "USA"
)

## solve the model with those parameters
solved <- solve_ode(pars)

## plot prevalence
vis_timeline(solved, "incidence", log = FALSE, freescales = TRUE)

## plot incidence
vis_timeline(solved, "incidence", log = FALSE, freescales = TRUE)



## DEFINE A LIST OF SCENARIOS FOR COMPARISON

## define a list of different scenarios
scen <- list(
  "No vaccination" = list(vax_rate = 0, vax_prioritised = FALSE),
  "Prioritised vaccination for elderly" = list(vax_rate = 0.002, vax_prioritised = TRUE),
  "Equal vaccination for all ages" = list(vax_rate = 0.002, vax_prioritised = FALSE)
) %>%
  map(~ solve_ode(get_pars(vax_rate = .x[[1]], vax_prioritised = .x[[2]])))

vis_comparison(scen, freescales = TRUE)
