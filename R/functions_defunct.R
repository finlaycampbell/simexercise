## define model parameters
get_pars_old <- function(R0=2.5,
                     generation_time = NULL, ## 8 days
                     doubling_time = 6.5, ## days
                     tau_E=5, ## age-stratified E duration
                     prop_E_trans=0.5, ## proportion of transmission from E
                     inf_E_over_I = 1, ## beta_E/beta_I,
                     frac_symp = 0.5, ## fraction symptomatic
                     symp_case_hosp_frac = c( ## fraction symptomatic hospitalized
                       0.0125,
                       rep(0.005,3),
                       rep(0.0125,6),
                       rep(0.0175, 3),
                       rep(0.16,3)
                     ),
                     rate_dschg = c( ## hospital discharge rate (1/duration of hospitalisations)
                       rep(1/8,4),
                       rep(1/9, 6),
                       rep(1/10, 6)
                     ),
                     max_day = 365,
                     tau_I = tau_E * inf_E_over_I * (1-prop_E_trans)/prop_E_trans,
                     ## get hospitalisation rate: this is essentially the odds
                     ## of entering hospital (ie for prob of 0.25 the odds are
                     ## 1/3) multiplied by the rate of leaving the I
                     ## compartment. So if you are 3x more likely to recover
                     ## than go to hospital, the hospitalisation rate is 1/3 of
                     ## the recovery rate (1/tau_I)
                     rate_hosp = (symp_case_hosp_frac/(1 - symp_case_hosp_frac))/tau_I,
                     scenario_type = "Baseline",
                     iso3 = "USA",
                     init_day = 1,
                     step_day = 1,
                     norm_100K = TRUE,
                     betas = TRUE
                     ) {

  ## form to list
  pars <- as.list(environment())

  if (!is.null(pars$generation_time) & is.null(pars$doubling_time)) {
    pars$doubling_time = with(pars, log(2) * generation_time/log(R0))
  }

  ## choose scenario parameters from potential scenarios
  pars$scenario <- scenarios[[scenario_type]]

  ## get age-group populations and age fractions
  pars$population = get_pop_frac(pop, iso3)$pop # Alt: get_pops_ACS(pars$GEOID)
  pars$age_frac = get_pop_frac(pop, iso3)$pop_frac # Proportion in age group

  ## get polymod matrix
  pars$poly <- get_poly(iso3)

  ## solve for transmission per contact using country specifi polymod matrix
  pars$p_trans =
    uniroot(function(x) trans_prob_to_doubling_time_old(x, pars) - pars$doubling_time, c(0, 1))$root

  ## normalize base population to 100K
  if (norm_100K) pars$population = 10^5 * pars$age_frac

  ## ## add
  ## if (betas) {
  ##   ## Default (baseline) matrices, modified for different scenarios:
  ##   pars$polymod = with(pars, t(t(poly$scale$all) * age_frac)) # Density-rescaling
  ##   ## Remember that beta_ij is for a 100% population of naive j
  ##   pars$beta_I = with(pars, p_trans * poly$scale$all)
  ##   pars$beta_E = with(pars, inf_E_over_I * beta_I)
  ## }

  ## ## Scaled polymods -- very important!
  ## polyscale.US = t(t(polymod.US)/prop.table(populations)) # m_ij/P_j
  ## polyscale.home.US = t(t(polymod.home.US)/prop.table(populations)) # m_ij/P_j
  ## polyscale.work.US = t(t(polymod.work.US)/prop.table(populations)) # m_ij/P_j
  ## polyscale.school.US = t(t(polymod.school.US)/prop.table(populations)) # m_ij/P_j
  ## polyscale.other.US = t(t(polymod.other.US)/prop.table(populations)) # m_ij/P_j
  ## polyscale.UNI = t(t(polymod.UNI)/prop.table(populations)) # m_ij/P_j
  ## polyscale.PK12 = t(t(polymod.PK12)/prop.table(populations)) # m_ij/P_j

  ## Symmetrize and population adjust polymod
  ## polymod.cuts = c(0:15 * 5, Inf)

  ## populations = get_pops16_ACS("US") # Use ACS data

  ## We don't seperate Uni and PK12 for now
  ## ## Separating out PK12 from university in school polymod (sweeping - but reasonable? - assumptions)
  ## ## Estimating polymod from polymod.school and population fractions in university
  ## ## Assume that p15to17 * school interactions 15-19 are at HS, remainder at university
  ## ## Assume all school interactions 20-24 are at universities
  ## ## Assume that only school interactions with 18-22 are fractionated out
  ## ## for the remainder of age groups as at university
  ## USUNI = get_popsUNI_ACS("US")
  ## p15to17 = (USUNI[1] - sum(populations[1:3]))/populations[4]
  ## p22to24 = (sum(populations[1:5]) - sum(USUNI[1:2]))/populations[5]
  ## polymod.UNI = c(0, 0, 0, 1-p15to17, p22to24, rep(0,11)) * populations * as.matrix(polymod.school.US)
  ## ## All 20-24 school interactions are at university (assumption):
  ## polymod.UNI[5,5] = populations[5] * polymod.school.US[5,5]
  ## ## Symmetrize number of 15-19 to 20-24 intercations:
  ## polymod.UNI[4,5] = (polymod.UNI[5,4] + polymod.UNI[4,5])/2
  ## polymod.UNI[5,4] = polymod.UNI[4,5]
  ## ## Symmetrize the rest:
  ## polymod.UNI[,4:5] = t(polymod.UNI[4:5,])
  ## ## Normalize to population
  ## polymod.UNI = polymod.UNI/populations
  ## polymod.PK12 = polymod.school.US - polymod.UNI

  ## calculate betas from pathogen characterstics, population structure and interventions
  pars <- add_betas_old(pars, mult = FALSE)

  return(pars)

}

## Solve the ODE
runODE_old <- function(pars) {

  ## Naive state: all susceptible
  state_naive <- cbind(
    S = pars$age_frac,
    E = 0*pars$age_frac,
    I = 0*pars$age_frac,
    R = 0*pars$age_frac
  ) %>%
    divide_by(sum(.)) %>%
    "row.names<-"(paste0("age_", seq_len(nrow(.))))

  ## Initial state may be supplied through pars
  if (is.null(pars$init_state)) {
    ## Initial state: seeded by default by a 35-39 year old (why not?)
    state <- state_naive
    midgp <- round(nrow(state)/2)
    ## 1 of 35-39 age group exposed
    state[midgp, c("S", "E")] <- state[midgp, "S"] * (c(1,0) + c(-1, 1)/pars$population[midgp])
  } else {
    state <- pars$init_state
  }

  ## Run ODE solver
  days <- seq(pars$init_day, pars$max_day, by=pars$step_day)
  solved <- ode(age_SEIR_old, y = c(state), times = days, parms=pars)

  ## Format output to age-stratified 3D array
  prevalence <- array(solved[,-1], dim=c(nrow(solved), dim(state_naive)))
  dimnames(prevalence) = list(
    time = paste0("day_", days),
    age = dimnames(state_naive)[[1]],
    state = dimnames(state_naive)[[2]]
  )

  ## calculate change in categories
  deltas <- vapply(
    seq_len(dim(prevalence)[1]),
    function(i) {
      if(i == 1) prevalence[i,,] - prevalence[i,,]
      else prevalence[i,,] - prevalence[i-1,,]
    },
    prevalence[1,,]
  ) %>%
    aperm(c(3, 1, 2)) %>%
    "dimnames<-"(dimnames(prevalence))

  ## incidence is the number of new entries into that category
  incidence <- array(
    c(
      ## S: we never have an entries into S
      rep(0, prod(dim(prevalence)[1:2])),
      ## E: infections of Susceptible from Exposed, and Susceptible from Infected
      t(apply(prevalence, 1, function(s) {
        (s[, "E"] %*% pars$beta_E) * s[, "S"] + (s[, "I"] %*% pars$beta_I) * s[, "S"]
      })),
      ## I: E multiplied by the rate of leaving E
      t(apply(prevalence, 1, function(s) s[,"E"]/pars$tau_E)),
      ## R: I multiplied by the rate of leaving I
      t(apply(prevalence, 1, function(s) s[,"I"]/pars$tau_I))
    ),
    dim = dim(prevalence),
    dimnames = dimnames(prevalence)
  )

  list(prevalence = prevalence, deltas = deltas, incidence = incidence, pars = pars)

}

## define ODEs for SEIR
age_SEIR_old <- function(Time, State, pars) {
  ## State is 4N-element vector (N-age-group X 4-SEIR)
  ## Each element is number of persons/total population (sum(State)=1)
  ## ODEs set up to work on N X 4 matrix
  nr = nrow(pars$beta_I) ## Gettin the number of age groups, N
  State = matrix(State, nrow=nr, dimnames=list(paste0("AGE", 1:nr), c("S", "E", "I", "R")))
  dState = State * 0
  dState[, "S"] = with(pars, - (State[, "E"] %*% beta_E) * State[, "S"] - (State[, "I"] %*% beta_I) * State[, "S"])
  dState[, "E"] = with(pars, - dState[, "S"] - State[, "E"]/tau_E)
  dState[, "I"] = with(pars, State[, "E"]/tau_E - State[, "I"]/tau_I)
  dState[, "R"] = with(pars, State[, "I"]/tau_I)
  list(dState)
}

## Converting from scenarios to betas using population structure
add_betas_old <- function(pars, mult=FALSE) {

  ## First combine the context specific pieces
  ## Adjusts for isolation duration being shorter than tau_I
  ## If !mult, more effective of layered intervention retained; else multiplicative effect
  ## pars$scenario is list with possible elements:
  ## home, school, work, other, PK12, UNI, isol, isol_eff, isol_dur, cocoon_eff

  ## ## split into PK12 and UNI if necessary
  ## Scenario$PK12 = max(Scenario$PK12, Scenario$school)
  ## Scenario$UNI = max(Scenario$UNI, Scenario$school)

  ## proportion of baseline contact rate under given social distancing scenario
  places <- c("home", "school", "work", "other")
  non_isol_coeffs <- unlist(pars$scenario[places])

  ## proportion of baseline contact rate under given isolation scenarios for
  ## home, PK12, uni, work, other (home is unaffected under isolation)
  ## this proportion is given by
  ## (proportion of sympomatic cases) x
  ## (proportion of symptomatic cases that are isolated) x
  ## (duration of isolation as proportion of symptomatic infectious period) x
  ## (1 - effectiveness of isolation)
  isol_coeffs <- c(
    1, ## for home
    rep(pars$frac_symp *
        pars$scenario$isol *
        max(pars$scenario$isol_dur/pars$tau_I) *
        (1 - pars$scenario$isol_eff), 3)
  )

  ## do social distancing and isolation have a combined effect?
  if (mult) isol_coeffs <- non_isol_coeffs * isol_coeffs
  else isol_coeffs <- pmin(non_isol_coeffs, isol_coeffs)

  ## calculation fraction of symptomatic cases that are isolated
  frac_isol <- pars$frac_symp * pars$scenario$isol * max(pars$scenario$isol_dur/pars$tau_I)

  ## Multiply non-isolation coefficients (i.e. proportion of contacts still
  ## occuring after a given intervention while *not* isolated) by baseline
  ## contact rates and sum to get total contacts in non-isolated state during a
  ## given intervention
  polyscale_non_isol <- map2(non_isol_coeffs, pars$poly$scale[places], \(x, y) x*y) %>%
    Reduce("+", .)

  ## Multiply isolation coefficients (i.e. proportion of contacts still occuring
  ## once isolated) by baseline contact rates and sum to get total contacts in
  ## isolated state
  polyscale_isol <- map2(isol_coeffs, pars$poly$scale[places], \(x, y) x*y) %>%
    Reduce("+", .)

  ## cocooning is a contact reduction specific to >65 year olds
  cocoon_coeffs = c(rep(1, 13), rep(1 - pars$scenario$cocoon_eff, 3))

  ## beta for Exposed cases takes:
  ## - contact rates unaffected by isolated but affected by social distancing
  ## - estimated transmission probability per contact
  ## - relative infectiousness of Exposed relative to Infected
  ## - cocooning effect on contact rates for elderly
  pars$beta_E = t(cocoon_coeffs * t(pars$inf_E_over_I * pars$p_trans * polyscale_non_isol))

  ## beta for Infected cases takes:
  ## - weighted mean of contact rates from isolated and non-isolated Infecteds
  ## - estimated transmission probability per contact
  ## - cocooning effect on contact rates for elderly
  pars$beta_I = t(cocoon_coeffs * t(pars$p_trans *
                                    (frac_isol * polyscale_isol + (1 - frac_isol) * polyscale_non_isol)))

  return(pars)

}

## convert transmission probability to doubling time
trans_prob_to_doubling_time_old <- function(p_trans, pars) {

  ## REM: beta_ij S_j I_i in naive population -> t(p_trans*poly.US) %*% I
  ## poly.US is polymod for the nation

  ## I think this is not the NGM (not multiplied by duration in categories?) but
  ## just the T/beta matrix?

  nr = length(pars$age_frac)
  b11 = t(p_trans*pars$inf_E_over_I*pars$poly$mod$all) - diag(1/pars$tau_E, nr)
  b12 = t(p_trans*pars$poly$mod$all)
  b21 = diag(1/pars$tau_E, nr)
  b22 = -diag(1/pars$tau_I, nr)
  b = rbind(cbind(b11, b12), cbind(b21, b22))
  res = log(2)/max(eigen(b)$values)
  ifelse(res<0, Inf, res)

}

## Post processing: hospitalized
get_hosp <- function(out_state, pars) {
  symp = cbind(Day=1:dim(out_state)[1], data.frame(t(apply(out_state[,, "I"], 1, function(rr) pars$frac_symp * rr))))
  hosp = symp[1,]
  hosp[, -1] = (1 - exp(-pars$rate_hosp*hosp$Day[1])) * hosp[, -1] # Zero on day 0
  hosp.symp = hosp # Hospitalized symptomatics; added 2020-03-31
  daily_adm = c()
  for (rr in 2:nrow(symp)) {
    adm = (symp[rr-1, -1] - hosp.symp[rr-1, -1]) * (1 - exp(-pars$rate_hosp*diff(symp$Day[rr-1:0])))
    not.dschg = hosp[rr-1, -1] * exp(-pars$rate_dschg*diff(symp$Day[rr-1:0]))
    not.dschg.symp = not.dschg * exp(-diff(symp$Day[rr-1:0])/pars$tau_I) # Not discharged and still symptomatic
    hosp.symp = rbind(hosp.symp, c(Day=symp$Day[rr], not.dschg.symp + adm))
    hosp = rbind(hosp, c(Day=symp$Day[rr], not.dschg + adm))
    daily_adm = rbind(daily_adm, adm)
  }
  list(hosp, daily_adm)
}

## code to batch process scenarios triggers and durations
allruns_scenario_duration_target <- function(geoid, social_distancing_scenarios, post_pause_scenario,
                                             closure_durations, targets, pars) {

  scenario.runs = list()
  ## Get populations (betas will be adjusted within runODE when polymod is modified)
  pars = pop_adjust_pars(pars, geoid, norm_100K=TRUE, betas = TRUE)
  for (target.i in 1:length(targets)) {
    cum_I_target = targets[target.i]
    allruns = list()
    ## Warnings ok:
    pars = within(pars, rm(init_state, init_day))

    ## Baseline run: no interventions
    pars$scenario = subset(social_distancing_scenarios, Scenario=="Baseline")

    ## Generate baserun -- duplicated across targets (redundant, but quick)
    baserun = runODE_old(pars)
    runs = colSums(baserun$inc_I)/sum(pars$population)
    runs_seniors = colSums(baserun$inc_I[14:16,])/sum(pars$population)
    run_log = data.frame(scenario="Baseline", duration=0)
    allruns[[nrow(run_log)]] = baserun$out_state
    runs_hosp = rowSums(get_hosp(baserun$out_state, pars)[[1]][, -1])
    cum_I = cumsum(colSums(baserun$inc_I)/sum(pars$population))
    ## Extract target date for intervention from baseline run
    target_dt = which(cum_I > cum_I_target)[1]

    ## Run interventions (social distancing scenarios)
    ##  First run for each scenario: indefinite intervention
    ##  Subsequent runs: social distancing stopped after different durations (longest first for easier computation)
    ##  Cocooning/isolation may be indefinitely imposed after termination of social distancing
    sc = setdiff(unique(social_distancing_scenarios$Scenario), "Baseline")
    for (sr in 1:length(sc)) {
      ## For each scenario, first an infinite duration run is generated
      pars$scenario = subset(social_distancing_scenarios, Scenario == sc[sr])
      pars$init_state = c(baserun$out_state[target_dt,,])
      pars$init_day = target_dt
      intrun = runODE_old(pars)
      int_inc_I = c(colSums(baserun$inc_I)[1:(target_dt - 1)], colSums(intrun$inc_I))/sum(pars$population)
      runs = rbind(runs, int_inc_I)
      runs_seniors = rbind(runs_seniors, c(colSums(baserun$inc_I[14:16,])[1:(target_dt - 1)],
                                           colSums(intrun$inc_I[14:16,]))/sum(pars$population))
      run_log = rbind(run_log, data.frame(scenario=social_distancing_scenarios$Scenario[sr], duration=Inf))
      allruns[[nrow(run_log)]] = allruns[[nrow(run_log)-1]]
      allruns[[nrow(run_log)]][target_dt:(dim(baserun$out_state)[1]),,] = intrun$out_state
      runs_hosp = rbind(runs_hosp, rowSums(get_hosp(allruns[[nrow(run_log)]], pars)[[1]][, -1]))

      ## Set parameters for post social distancing phase
      pars$scenario = post_pause_scenario

      ## Running in decreasing order of closure duration (allows each run to be based off segment of the previous one)
      closure_durations = sort(closure_durations, decreasing=TRUE)
      for (int_duration in closure_durations) { # Longest duration first, as log borrows from previous run
        browser()
        pars$init_state = c(intrun$out_state[int_duration + 1,,])
        end_dt = target_dt + int_duration # end of intervention
        pars$init_day = end_dt
        int_end_run = runODE_old(pars)
        end_inc_I = c(tail(runs[,1:(end_dt-1)], 1), colSums(int_end_run$inc_I)/sum(pars$population))
        runs = rbind(runs, end_inc_I)
        runs_seniors = rbind(runs_seniors,c(tail(runs_seniors[,1:(end_dt-1)], 1), colSums(int_end_run$inc_I[14:16,])/sum(pars$population)))
        run_log = rbind(run_log, data.frame(scenario=sc[sr], duration=int_duration))
        allruns[[nrow(run_log)]] = allruns[[nrow(run_log)-1]]
        allruns[[nrow(run_log)]][end_dt:(dim(baserun$out_state)[1]),,] = int_end_run$out_state
        runs_hosp = rbind(runs_hosp, rowSums(get_hosp(allruns[[nrow(run_log)]], pars)[[1]][, -1]))
      }
    }

    ## Delete some parameters (not really necessary):
    pars = within(pars, rm(init_state, init_day))
    ## Row name warning ok:
    peaks = cbind(run_log, cum_I_target=cum_I_target, target_dt=target_dt, end_dt=target_dt+ run_log$duration,
                  peak_date = apply(runs, 1, which.max), peak_incidence_100K = 10^5*apply(runs, 1, max),
                  peak_date_seniors = apply(runs_seniors, 1, which.max), peak_incidence_seniors_100K = 10^5*apply(runs_seniors, 1, max),
                  peak_date_hosp = apply(runs_hosp, 1, which.max), peak_hosp_100K = 10^5*apply(runs_hosp, 1, max),
                  post_trigger_week_cases = 10^5 * rowSums(runs[, target_dt[1]+1:7]))
    scenario.runs[[target.i]] = list(peaks=peaks, runs=runs, runs_seniors=runs_seniors, runs_hosp=runs_hosp, allruns=allruns)
  }
  scenario.runs
}
