setwd("c:/users/shubho/covid")
options(stringsAsFactors=FALSE)
library(rio)
library(deSolve)

# Country specific social contact rate data (suppl to Prem et al)
polymod = import("contact_matrices_152_countries/MUestimates_all_locations_2.xlsx",
  which="United States of America", col_names=FALSE)
polymod.home = import("contact_matrices_152_countries/MUestimates_home_2.xlsx",
  which="United States of America", col_names=FALSE)
polymod.school = import("contact_matrices_152_countries/MUestimates_school_2.xlsx",
  which="United States of America", col_names=FALSE)
polymod.work = import("contact_matrices_152_countries/MUestimates_work_2.xlsx",
  which="United States of America", col_names=FALSE)
polymod.other = import("contact_matrices_152_countries/MUestimates_other_locations_2.xlsx",
  which="United States of America", col_names=FALSE)

get_pops16_ACS = function(geoid) { # Function to get population in 16 age categories
 rep(10^8, 16)
}

get_popsUNI_ACS = function(geoid = "13089") { # High-school, university and older populations
 # 3-category population: 0-17, 18-21, 22+
  pops = c(4E8, 1E8, 11E8)
  names(pops) = c("HS", "UNI", "OLD")
  pops
}

trans_prob_to_doubling_time = function(trans_per_cont, Pars) {
  # REM: beta_ij S_j I_i in naive population -> t(trans_per_cont*poly.US) %*% I
  # poly.US is polymod for the nation
  nr = length(Pars$age_frac)
  b11 = t(trans_per_cont*Pars$inf_E_over_I*poly.US) - diag(1/Pars$tau_E, nr)
  b12 = t(trans_per_cont*poly.US)
  b21 = diag(1/Pars$tau_E, nr)
  b22 = -diag(1/Pars$tau_I, nr)
  b = rbind(cbind(b11, b12), cbind(b21, b22))
  res = log(2)/max(eigen(b)$values)
  ifelse(res<0, Inf, res)
}

# Doubling time and principal eigenvector from next generation matrix
next_gen_analysis = function(Pars, verbose=FALSE) {
  with(Pars, {
    nr = nrow(beta_I)
    b11 = age_frac*t(beta_E) - diag(1/tau_E, nr)
    b12 = age_frac*t(beta_I)
    b21 = diag(1/tau_E, nr)
    b22 = -diag(1/tau_I, nr)
    b = rbind(cbind(b11, b12), cbind(b21, b22))
    ev = eigen(b)
    mval = max(Re(ev$values))
    mvec = ev$vectors[, which.max(Re(ev$values))]
    check.sign = (sum(abs(sign(Re(mvec)))) == abs(sum(sign(Re(mvec))))) # Check if max vector is all of same sign or zero
    if (!check.sign & verbose) print("Anomalous eigenvector of largest eigenvalue detected")
    list(doubling_time=log(2)/max(0, Re(eigen(b)$values)), eval=ev$values[which.max(Re(ev$values))], evec=mvec)
  })
}

pop_adjust_pars = function(Pars, geoid="US", norm_100K=FALSE, betas=FALSE) {
  # Boiler plate updates
  Pars$max_day = 365 # Simulation period
  Pars$tau_I = with(Pars,  tau_E * inf_E_over_I * (1-prop_E_trans)/prop_E_trans)
  Pars$rate_hosp = with(Pars, (symp_case_hosp_frac/(1 - symp_case_hosp_frac))/tau_I) # get_hosp applies rate to frac_sym, 1/(1-.) correction 2020-03-31
  if (!is.null(Pars$generation_time) & is.null(Pars$doubling_time)) {
    Pars$doubling_time = with(Pars, log(2) * generation_time/log(R0))
  }

  Pars$GEOID = geoid
  Pars$population = get_pops16_ACS(Pars$GEOID) # Alt: get_pops_ACS(Pars$GEOID)
  Pars$age_frac = prop.table(Pars$population) # Proportion in age group

  poly.US <<- polymod.US
  Pars$prob_trans_per_contact = uniroot(function(x) trans_prob_to_doubling_time(x, Pars)-Pars$doubling_time, c(0, 1))$root


  if (norm_100K) Pars$population = 10^5 * Pars$age_frac # Normalize base population to 100K
  
  if (betas) {
    # Default (baseline) matrices, modified for different scenarios:
    Pars$polymod = with(Pars, t(t(polyscale.US) * age_frac)) # Density-rescaling
    # Remember that beta_ij is for a 100% population of naive j
    Pars$beta_I = with(Pars, prob_trans_per_contact * polyscale.US)
    Pars$beta_E = with(Pars, inf_E_over_I * beta_I)
  }

  Pars
}

# Post processing: hospitalized
get_hosp = function(out_state, Pars) {
    symp = cbind(Day=1:dim(out_state)[1], data.frame(t(apply(out_state[,, "I"], 1, function(rr) Pars$frac_symp * rr))))
    hosp = symp[1,]
    hosp[, -1] = (1 - exp(-Pars$rate_hosp*hosp$Day[1])) * hosp[, -1] # Zero on day 0
    hosp.symp = hosp # Hospitalized symptomatics; added 2020-03-31
    daily_adm = c()
    for (rr in 2:nrow(symp)) {
      adm = (symp[rr-1, -1] - hosp.symp[rr-1, -1]) * (1 - exp(-Pars$rate_hosp*diff(symp$Day[rr-1:0])))
      not.dschg = hosp[rr-1, -1] * exp(-Pars$rate_dschg*diff(symp$Day[rr-1:0]))
      not.dschg.symp = not.dschg * exp(-diff(symp$Day[rr-1:0])/Pars$tau_I) # Not discharged and still symptomatic
      hosp.symp = rbind(hosp.symp, c(Day=symp$Day[rr], not.dschg.symp + adm))
      hosp = rbind(hosp, c(Day=symp$Day[rr], not.dschg + adm))
      daily_adm = rbind(daily_adm, adm)
    }
    list(hosp, daily_adm)
}



## ----defining the main ODE solving wrappers, warning=FALSE, message=FALSE, echo=TRUE------------------------------------------------------------------

age_SEIR = function(Time, State, Pars) {
  # State is 4N-element vector (N-age-group X 4-SEIR)
  # Each element is number of persons/total population (sum(State)=1)
  # ODEs set up to work on N X 4 matrix
  nr = nrow(Pars$beta_I) # Gettin the number of age groups, N
  State = matrix(State, nrow=nr, dimnames=list(paste0("AGE", 1:nr), c("S", "E", "I", "R"))) 
  dState = State * 0
  dState[, "S"] = with(Pars, - (State[, "E"] %*% beta_E) * State[, "S"] - (State[, "I"] %*% beta_I) * State[, "S"])
  dState[, "E"] = with(Pars, - dState[, "S"] - State[, "E"]/tau_E)
  dState[, "I"] = with(Pars, State[, "E"]/tau_E - State[, "I"]/tau_I)
  dState[, "R"] = with(Pars, State[, "I"]/tau_I)
  list(dState)
}

### Converting from scenarios to betas
scenario2beta = function(Pars, mult=FALSE) {
  ## First combine the context specific pieces
  ## Adjusts for isolation duration being shorter than tau_I
  ## If !mult, more effective of layered intervention retained; else multiplicative effect
  ## Pars$scenario is list with possible elements:
  ## home, school, work, other, PK12, UNI, isol, isol_eff, isol_dur, cocoon_eff
  Scenario = as.list(Pars$scenario)
  Scenario$PK12 = max(Scenario$PK12, Scenario$school)
  Scenario$UNI = max(Scenario$UNI, Scenario$school)
  with(as.list(Scenario), {
    non_isol_coeffs = c(home, PK12, UNI, work, other)
    isol_coeffs = c(1, rep(Pars$frac_symp * isol * max(isol_dur/Pars$tau_I) * (1 - isol_eff), 4))
    if (mult) {isol_coeffs = non_isol_coeffs * isol_coeffs} else {isol_coeffs = pmin(non_isol_coeffs, isol_coeffs)}
    frac_isol = Pars$frac_symp * isol * max(isol_dur/Pars$tau_I)
    polyscale_non_isol = non_isol_coeffs[1] * polyscale.home.US + non_isol_coeffs[2] * polyscale.PK12 +
      non_isol_coeffs[3] * polyscale.UNI + non_isol_coeffs[4] * polyscale.work.US + non_isol_coeffs[5] * polyscale.other.US
    polyscale_isol = isol_coeffs[1] * polyscale.home.US + isol_coeffs[2] * polyscale.PK12 +
      isol_coeffs[3] * polyscale.UNI + isol_coeffs[4] * polyscale.work.US + isol_coeffs[5] * polyscale.other.US
    cocoon_coeffs = c(rep(1, 13), rep(1 - cocoon_eff, 3))
    Pars$beta_E = t(cocoon_coeffs * t(Pars$inf_E_over_I * Pars$prob_trans_per_contact * polyscale_non_isol))
    Pars$beta_I = t(cocoon_coeffs * t(Pars$prob_trans_per_contact * (frac_isol * polyscale_isol + (1 - frac_isol) * polyscale_non_isol)))
    list(Pars=Pars, non_isol_coeffs=non_isol_coeffs, isol_coeffs=isol_coeffs, frac_isol=frac_isol)
  })
}

runODE = function(Pars) {
  # Naive state: all susceptible
  state_naive = cbind(S=Pars$age_frac, E=0*Pars$age_frac, I=0*Pars$age_frac, R=0*Pars$age_frac)
  state_naive = state_naive/sum(state_naive) # For good measure

  # Initial state may be supplied through Pars
  if (is.null(Pars$init_state)) {
    # Initial state: seeded by default by a 35-39 year old (why not?)
    state = state_naive
    midgp = round(length(Pars$age_frac)/2)
    state[midgp, c("S", "E")] = state[midgp, "S"] * (c(1,0) + c(-1, 1)/Pars$population[midgp]) # 1 of 35-39 age group exposed
  } else {
    state = Pars$init_state
  }

  Pars = scenario2beta(Pars)$Pars

  if (is.null(Pars$init_day)) Pars$init_day = 1
  if (is.null(Pars$step_day)) Pars$step_day = 1 # step size, for output

  # Run ODE solver
  days = seq(Pars$init_day, Pars$max_day, by=Pars$step_day)
  out = ode(age_SEIR, y = c(state), times = days, parms=Pars)
  
  # Format output to age-stratified 3D array
  out_state = array(out[,-1], dim=c(nrow(out), dim(state_naive)))
  dimnames(out_state) = list(time=paste("Day", out[,1]), age=dimnames(state_naive)[[1]], state=dimnames(state_naive)[[2]])

  # Post-processing: incident cases (daily count) of I and E
  inc_I = apply(out_state, 1, function(s) s[,"E"]/Pars$tau_E) * sum(Pars$population)
  inc_E = apply(out_state, 1, function(s) {
    (s[, "E"] %*% Pars$beta_E) * s[, "S"] + (s[, "I"] %*% Pars$beta_I) * s[, "S"]
    }) * sum(Pars$population)
  peak_incidence = max(colSums(matrix(inc_I, nrow=dim(out_state)[2]))) # Ensuring matrix for single age-group models
  peak_incidence_date = days[which.max(colSums(matrix(inc_I, nrow=dim(out_state)[2])))]

  list(out_state=out_state, inc_E=inc_E, inc_I=inc_I, peak_incidence=peak_incidence, peak_incidence_date=peak_incidence_date, 
    Pars=Pars)
}



## ----code to batch process scenarios triggers and durations, warning=FALSE, message=FALSE-------------------------------------------------------------
allruns_scenario_duration_target = function(geoid, social_distancing_scenarios, post_pause_scenario, closure_durations, targets, Pars) {
  scenario.runs = list()
  # Get populations (betas will be adjusted within runODE when polymod is modified)
  Pars = pop_adjust_pars(Pars, geoid, norm_100K=TRUE)
  for (target.i in 1:length(targets)) {
    cum_I_target = targets[target.i]
    allruns = list()
    # Warnings ok:
    Pars = within(Pars, rm(init_state, init_day))
    
    ## Baseline run: no interventions
    Pars$scenario = subset(social_distancing_scenarios, Scenario=="Baseline")

    # Generate baserun -- duplicated across targets (redundant, but quick)
    baserun = runODE(Pars)
    runs = colSums(baserun$inc_I)/sum(Pars$population)
    runs_seniors = colSums(baserun$inc_I[14:16,])/sum(Pars$population)
    run_log = data.frame(scenario="Baseline", duration=0)
    allruns[[nrow(run_log)]] = baserun$out_state
    runs_hosp = rowSums(get_hosp(baserun$out_state, Pars)[[1]][, -1])
    cum_I = cumsum(colSums(baserun$inc_I)/sum(Pars$population))
    # Extract target date for intervention from baseline run
    target_dt = which(cum_I > cum_I_target)[1]

    ## Run interventions (social distancing scenarios)
    ##  First run for each scenario: indefinite intervention
    ##  Subsequent runs: social distancing stopped after different durations (longest first for easier computation)
    ##  Cocooning/isolation may be indefinitely imposed after termination of social distancing 
    sc = setdiff(unique(social_distancing_scenarios$Scenario), "Baseline")
    for (sr in 1:length(sc)) { 
      # For each scenario, first an infinite duration run is generated
      Pars$scenario = subset(social_distancing_scenarios, Scenario == sc[sr])
      Pars$init_state = c(baserun$out_state[target_dt,,])
      Pars$init_day = target_dt
      intrun = runODE(Pars)
      int_inc_I = c(colSums(baserun$inc_I)[1:(target_dt - 1)], colSums(intrun$inc_I))/sum(Pars$population)
      runs = rbind(runs, int_inc_I)
      runs_seniors = rbind(runs_seniors, c(colSums(baserun$inc_I[14:16,])[1:(target_dt - 1)],
        colSums(intrun$inc_I[14:16,]))/sum(Pars$population))
      run_log = rbind(run_log, data.frame(scenario=social_distancing_scenarios$Scenario[sr], duration=Inf))
      allruns[[nrow(run_log)]] = allruns[[nrow(run_log)-1]]
      allruns[[nrow(run_log)]][target_dt:(dim(baserun$out_state)[1]),,] = intrun$out_state
      runs_hosp = rbind(runs_hosp, rowSums(get_hosp(allruns[[nrow(run_log)]], Pars)[[1]][, -1]))

      # Set parameters for post social distancing phase
      Pars$scenario = post_pause_scenario
 
      # Running in decreasing order of closure duration (allows each run to be based off segment of the previous one)
      closure_durations = sort(closure_durations, decreasing=TRUE)
      for (int_duration in closure_durations) { # Longest duration first, as log borrows from previous run
        Pars$init_state = c(intrun$out_state[int_duration + 1,,])
        end_dt = target_dt + int_duration # end of intervention
        Pars$init_day = end_dt
        int_end_run = runODE(Pars)
        end_inc_I = c(tail(runs[,1:(end_dt-1)], 1), colSums(int_end_run$inc_I)/sum(Pars$population))
        runs = rbind(runs, end_inc_I)
        runs_seniors = rbind(runs_seniors,c(tail(runs_seniors[,1:(end_dt-1)], 1), colSums(int_end_run$inc_I[14:16,])/sum(Pars$population)))
        run_log = rbind(run_log, data.frame(scenario=sc[sr], duration=int_duration))
        allruns[[nrow(run_log)]] = allruns[[nrow(run_log)-1]]
        allruns[[nrow(run_log)]][end_dt:(dim(baserun$out_state)[1]),,] = int_end_run$out_state
        runs_hosp = rbind(runs_hosp, rowSums(get_hosp(allruns[[nrow(run_log)]], Pars)[[1]][, -1]))
      }
  }
  
  # Delete some parameters (not really necessary):
  Pars = within(Pars, rm(init_state, init_day))
  # Row name warning ok:
  peaks = cbind(run_log, cum_I_target=cum_I_target, target_dt=target_dt, end_dt=target_dt+ run_log$duration,
    peak_date = apply(runs, 1, which.max), peak_incidence_100K = 10^5*apply(runs, 1, max),
    peak_date_seniors = apply(runs_seniors, 1, which.max), peak_incidence_seniors_100K = 10^5*apply(runs_seniors, 1, max),
    peak_date_hosp = apply(runs_hosp, 1, which.max), peak_hosp_100K = 10^5*apply(runs_hosp, 1, max),
    post_trigger_week_cases = 10^5 * rowSums(runs[, target_dt[1]+1:7]))
  scenario.runs[[target.i]] = list(peaks=peaks, runs=runs, runs_seniors=runs_seniors, runs_hosp=runs_hosp, allruns=allruns)
  }
  scenario.runs
}

## End helper functions



## ----Symmetrize and population adjust polymod, warning=FALSE, message=FALSE---------------------------------------------------------------------------
#######################################################
polymod.cuts = c(0:15 * 5, Inf)
US16 = get_pops16_ACS("US") # Use ACS data
# Regularize the darned polymod matrix for the national population distribution
# Arithmetic mean, to ensure additivity and consistency with flu-models approach
polymod.US = ((US16 * polymod) + t(US16 * polymod))/(2 * US16)
polymod.school.US = ((US16 * polymod.school) + t(US16 * polymod.school))/(2 * US16)
polymod.work.US = ((US16 * polymod.work) + t(US16 * polymod.work))/(2 * US16)
polymod.home.US = ((US16 * polymod.home) + t(US16 * polymod.home))/(2 * US16)
polymod.other.US = ((US16 * polymod.other) + t(US16 * polymod.other))/(2 * US16)

# Separating out PK12 from university in school polymod (sweeping - but reasonable? - assumptions)
# Estimating polymod from polymod.school and population fractions in university
# Assume that p15to17 * school interactions 15-19 are at HS, remainder at university
# Assume all school interactions 20-24 are at universities
# Assume that only school interactions with 18-22 are fractionated out for the remainder of age groups as at university
USUNI = get_popsUNI_ACS("US")
p15to17 = (USUNI[1] - sum(US16[1:3]))/US16[4]
p22to24 = (sum(US16[1:5]) - sum(USUNI[1:2]))/US16[5]
polymod.UNI = c(0, 0, 0, 1-p15to17, p22to24, rep(0,11)) * US16 * as.matrix(polymod.school.US)
# All 20-24 school interactions are at university (assumption):
polymod.UNI[5,5] = US16[5] * polymod.school.US[5,5]
# Symmetrize number of 15-19 to 20-24 intercations:
polymod.UNI[4,5] = (polymod.UNI[5,4] + polymod.UNI[4,5])/2
polymod.UNI[5,4] = polymod.UNI[4,5]
# Symmetrize the rest:
polymod.UNI[,4:5] = t(polymod.UNI[4:5,])
# Normalize to population
polymod.UNI = polymod.UNI/US16
polymod.PK12 = polymod.school.US - polymod.UNI

## Scaled polymods -- very important!
polyscale.US = t(t(polymod.US)/prop.table(US16)) # m_ij/P_j
polyscale.home.US = t(t(polymod.home.US)/prop.table(US16)) # m_ij/P_j
polyscale.work.US = t(t(polymod.work.US)/prop.table(US16)) # m_ij/P_j
polyscale.school.US = t(t(polymod.school.US)/prop.table(US16)) # m_ij/P_j
polyscale.other.US = t(t(polymod.other.US)/prop.table(US16)) # m_ij/P_j
polyscale.UNI = t(t(polymod.UNI)/prop.table(US16)) # m_ij/P_j
polyscale.PK12 = t(t(polymod.PK12)/prop.table(US16)) # m_ij/P_j



## ----set up parameters, warning=FALSE, message=FALSE, echo=TRUE---------------------------------------------------------------------------------------
#######################################################
# Set up parameters
# Updated 2020-04-02
if (!("Proposed Outbreak Scenarios COVID 19 current.rds" %in% list.files())) {
pars_scenarios = list(
 Best = list(R0=2.5, 
            # generation_time = 7, # days;
            doubling_time = 5.5, # days
            tau_E=5, # age-stratified E duration
            prop_E_trans=0.35, # proportion of transmission from E
            inf_E_over_I = 1, # beta_E/beta_I,
            frac_symp = 0.65,
            symp_case_hosp_frac = c(0.045, rep(0.006,3), rep(0.035,6), rep(0.07, 3), rep(0.1,3)), # fraction symptomatic hospitalized 
            rate_dschg = c(rep(1/3.3,10), rep(1/4, 3), rep(1/5.1,3)) # hospital discharge rate
            ),
 S1 = list(R0=2, 
            # generation_time = 7, # days
            doubling_time = 7, # days
            tau_E=5, # age-stratified E duration
            prop_E_trans=0.20, # proportion of transmission from E
            inf_E_over_I = 0.5, # beta_E/beta_I,
            frac_symp = 0.8,
            symp_case_hosp_frac = c(0.005, rep(0.001,3), rep(0.02,6), rep(0.045, 3), rep(0.045,3)), # fraction symptomatic hospitalized 
            rate_dschg = c(rep(1/3.3,10), rep(1/4, 3), rep(1/5.1,3)) # hospital discharge rate
            ),
 S2 = list(R0=2, 
            # generation_time = 8, # days
            doubling_time = 8, # days
            tau_E=5, # age-stratified E duration
            prop_E_trans=0.40, # proportion of transmission from E
            inf_E_over_I = 1, # beta_E/beta_I,
            frac_symp = 0.5,
            symp_case_hosp_frac = c(0.005, rep(0.001,3), rep(0.02,6), rep(0.045, 3), rep(0.045,3)), # fraction symptomatic hospitalized 
            rate_dschg = c(rep(1/3.3,10), rep(1/4, 3), rep(1/5.1,3)) # hospital discharge rate
            ),
 S3 = list(R0=3, 
            # generation_time = 8, # days
            doubling_time = 5, # days
            tau_E=5, # age-stratified E duration
            prop_E_trans=0.20, # proportion of transmission from E
            inf_E_over_I = 0.5, # beta_E/beta_I,
            frac_symp = 0.8,
            symp_case_hosp_frac = c(0.12, rep(0.015,3), rep(0.065,6), rep(0.1, 3), rep(0.17,3)), # fraction symptomatic hospitalized 
            rate_dschg = c(rep(1/3.3,10), rep(1/4, 3), rep(1/5.1,3)) # hospital discharge rate
            ),
 S4 = list(R0=3, 
            # generation_time = 8, # days
            doubling_time = 5, # days
            tau_E=5, # age-stratified E duration
            prop_E_trans=0.40, # proportion of transmission from E
            inf_E_over_I = 1, # beta_E/beta_I,
            frac_symp = 0.5,
            symp_case_hosp_frac = c(0.12, rep(0.015,3), rep(0.065,6), rep(0.1, 3), rep(0.17,3)), # fraction symptomatic hospitalized 
            rate_dschg = c(rep(1/3.3,10), rep(1/4, 3), rep(1/5.1,3)) # hospital discharge rate
            )
)
saveRDS(pars_scenarios, "Proposed Outbreak Scenarios COVID 19 current.rds")
} else { pars_scenarios = readRDS("Proposed Outbreak Scenarios COVID 19 current.rds")}

# Set population adjusted parameters
pars = pars_scenarios$Best
pars = pop_adjust_pars(pars, "US")



## ----Define standard scenarios to compare, warning=FALSE, message=FALSE-------------------------------------------------------------------------------
##########################################################################################
# Some "standard" social distancing scenarios
social_distancing_scenarios = data.frame(Scenario="Baseline", home=1, school=1, work=1, other=1, cocoon_eff=0, isol=0, isol_eff=0)
social_distancing_scenarios = rbind(social_distancing_scenarios, 
  data.frame(Scenario="School closure", home=1, school=0.2, work=1, other=1, cocoon_eff=0, isol=0, isol_eff=0),
  data.frame(Scenario="Telework 50%", home=1, school=1, work=0.5, other=1, cocoon_eff=0, isol=0, isol_eff=0),
  data.frame(Scenario="Home isolation", home=1, school=1, work=1, other=1, cocoon_eff=0, isol=0.5, isol_eff=0.75),
  data.frame(Scenario="School closure, telework 50%", home=1, school=0.2, work=0.5, other=1, cocoon_eff=0, isol=0, isol_eff=0),
  data.frame(Scenario="Lockdown or pause, home isolation", home=1, school=0.2, work=0.5, other=0.5, cocoon_eff=0, isol=0.5, isol_eff=0.75),
  data.frame(Scenario="Cocooning", home=1, school=1, work=1, other=1, cocoon_eff=0.5, isol=0, isol_eff=0),
  data.frame(Scenario="Community distancing, home isolation", home=1, school=1, work=0.5, other=0.5, cocoon_eff=0, isol=0.5, isol_eff=0.75),
  data.frame(Scenario="Community distancing, home isolation, cocooning", home=1, school=1, work=0.5, other=0.5, cocoon_eff=0.5, isol=0.5, isol_eff=0.75),
  data.frame(Scenario="General distancing, home isolation", home=1, school=0.8, work=0.8, other=0.8, cocoon_eff=0, isol=0.5, isol_eff=0.75),      
  data.frame(Scenario="Curfew, home isolation", home=1, school=0.2, work=0.2, other=0.2, cocoon_eff=0, isol=0.5, isol_eff=0.75)
)
social_distancing_scenarios$isol_dur = with(social_distancing_scenarios, 7 * (isol>0)) # Not necessary, bu less confusing



## ----large bundle of scenarious by trigger and duration, warning=FALSE, message=FALSE-----------------------------------------------------------------
##########################################################################
# Detailed analysis of scenarios with trigger and duration
# Try school closure, etc, triggered at cum_I_target cumulative I per 100K, for various lengths of time
## Big bundle of runs in scenario.runs for easy comparison
geoid = "US"
targets = c(100, 1000, 10000)/10^5
closure_durations = c(6, 8, 12, 16) * 7 # 
all_sc = unique(social_distancing_scenarios$Scenario)
for (pps in all_sc) {
  post_pause_scenario = subset(social_distancing_scenarios, Scenario==pps)
  pars_sc = "Best"

  # Warnings expected:
  myruns = allruns_scenario_duration_target(geoid, social_distancing_scenarios, post_pause_scenario, closure_durations, targets, pars_scenarios[[pars_sc]])
  out_data = list(geoid=geoid, social_distancing_scenarios=social_distancing_scenarios, post_pause_scenario=post_pause_scenario, closure_durations=closure_durations,
    targets=targets, pars=pars_scenarios[[pars_sc]], scenario.runs=myruns)
  saveRDS(out_data, paste0("out_data_", geoid, "_", pars_sc, "_", gsub("[,%+]", "", post_pause_scenario$Scenario[1]), "_",Sys.Date(), ".rds"))
}


