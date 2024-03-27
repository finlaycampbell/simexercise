## define model parameters
get_pars <- function(R0 = 2.5,
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
                     iso3 = "USA",
                     init_day = 1,
                     step_day = 1,
                     norm_100K = TRUE
                     ) {

  ## form to list
  pars <- as.list(environment())

  ## choose scenario parameters from potential scenarios
  pars$scenario <- scenarios[[scenario_type]]

  ## get age-group populations and age fractions
  pars$population = cdat[[iso3]]$pop$count
  pars$age_frac = cdat[[iso3]]$pop$prop

  ## get polymod contact matrix
  pars$poly <- list(mod = cdat[[iso3]]$mod, scale = cdat[[iso3]]$scale)

  if(generation_time < incubation_period/2)
    stop("Generation time not possible with incubation period provided.")

  ## calculate mortality if you would go to hospital but can't
  pars$unhosp_mortality <- hosp_mortality/(1-hosp_protection)

  ## solve for infectious period in community to satisfy generation time
  ## (incubation period and hospital period are fixed)
  pars$comm_duration <- uniroot(
    function(comm_duration) {
      pars$comm_duration <- comm_duration
      ## infection prob does not affect gen time so pick arbitrary one
      get_generation_time(0.01, add_rates(pars)) - generation_time
    },
    interval = c(0.01, generation_time*5)
  )$root

  ## add rates with correct comm_duration
  pars %<>% add_rates()

  ## calibrate p_trans from R0 or doubling time
  pars$p_trans <- uniroot(\(x) get_R0(x, pars) - R0, c(0, 1))$root
  pars$doubling_time <- get_doubling_time(pars$p_trans, pars)

  ## normalize base population to 100K
  if (norm_100K) pars$population = 10^5 * pars$age_frac

  ## calculate betas from pathogen, population structure and interventions
  pars %<>% add_betas(mult = FALSE)

  return(pars)

}

## Solve the ODE
solve_ode <- function(pars) {

  ## Naive state: all susceptible
  state_naive <- cbind(
    S_u = pars$age_frac,
    E_u = 0*pars$age_frac,
    C_u = 0*pars$age_frac,
    H_u = 0*pars$age_frac,
    R_u = 0*pars$age_frac,
    D_u = 0*pars$age_frac,
    S_v = 0*pars$age_frac,
    E_v = 0*pars$age_frac,
    C_v = 0*pars$age_frac,
    H_v = 0*pars$age_frac,
    R_v = 0*pars$age_frac,
    D_v = 0*pars$age_frac
  ) %>%
    divide_by(sum(.)) %>%
    "row.names<-"(paste0("age_", seq_len(nrow(.))))

  ## Initial state may be supplied through pars
  if (is.null(pars$init_state)) {
    ## Initial state: seeded by default by a 35-39 year old (why not?)
    state <- state_naive
    midgp <- round(nrow(state)/2)
    ## 1 of 35-39 age group exposed
    state[midgp, c("S_u", "E_u")] <-
      state[midgp, "S_u"] * (c(1,0) + c(-1, 1)/pars$population[midgp])
  } else {
    state <- pars$init_state
  }

  ## Run ODE solver
  days <- seq(pars$init_day, pars$max_day, by=pars$step_day)
  solved <- ode(age_SEIR, y = c(state), times = days, parms=pars)

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
      ## E_u: infections of Susceptible from Exposed, and Susceptible from Infected
      t(apply(prevalence, 1, \(state) {
        state[, "S_u"] * (
          state[, "E_u"] %*% pars$beta_E  +
          state[, "C_u"] %*% pars$beta_I_c +
          state[, "H_u"] %*% pars$beta_I_h +
          (1-pars$vax_infectious) * state[, "E_v"] %*% pars$beta_E +
          (1-pars$vax_infectious) * state[, "C_v"] %*% pars$beta_I_c +
          (1-pars$vax_infectious) * state[, "H_v"] %*% pars$beta_I_h
        )
      })),
      ## C_u: E_u multiplied by the rate of leaving E_u into cu
      t(apply(prevalence, 1, \(state) state[,"E_u"]*pars$lambda_cu)),
      ## H_u: E_u multiplied by the rate of leaving E_u into hu
      t(apply(prevalence, 1, \(state) state[,"E_u"]*pars$lambda_hu)),
      ## R_u: I_c * community recovery rate + I_h * hosp recovery rate (unvax)
      t(apply(prevalence, 1, \(state) {
        state[,"C_u"]*pars$sigma_cu + state[,"H_u"]*pars$sigma_hu
      })),
      ## D_u: I_c * community mortality rate + I_h * hosp mortality rate (unvax)
      t(apply(prevalence, 1, \(state) {
        state[,"C_u"]*pars$mu_cu + state[,"H_u"]*pars$mu_hu
      })),
      ## S_v: entry from S
      t(apply(prevalence, 1, \(state) state[,"S_u"]*pars$vax_rate)),
      ## E_v: V being infected by E_u, C_u and H_u, E_v, C_v, H_v
      ## these are protected from infection by vax_infection
      ## infectiousness reduction of vaccinated incorporated here
      t(apply(prevalence, 1, \(state) {
        state[, "S_v"] * (1-pars$vax_infection) * (
          state[, "E_u"] %*% pars$beta_E  +
          state[, "C_u"] %*% pars$beta_I_c +
          state[, "H_u"] %*% pars$beta_I_h +
          (1-pars$vax_infectious) * state[, "E_v"] %*% pars$beta_E +
          (1-pars$vax_infectious) * state[, "C_v"] %*% pars$beta_I_c +
          (1-pars$vax_infectious) * state[, "H_v"] %*% pars$beta_I_h
        )
      })),
      ## C_v: E_v multiplied by the rate of leaving E_v into cv
      t(apply(prevalence, 1, \(state) state[,"E_v"]*pars$lambda_cv)),
      ## H_v: E_v multiplied by the rate of leaving E_v into hv
      t(apply(prevalence, 1, \(state) state[,"E_v"]*pars$lambda_hv)),
      ## R_v: I_c * community recovery rate + I_h * hosp recovery rate (for vax)
      t(apply(prevalence, 1, \(state) {
        state[,"C_v"]*pars$sigma_cv + state[,"H_v"]*pars$sigma_hv
      })),
      ## D_v: I_c * community mortality rate + I_h * hosp mortality rate (for vax)
      t(apply(prevalence, 1, \(state) {
        state[,"C_v"]*pars$mu_cv + state[,"H_v"]*pars$mu_hv
      }))
    ),
    dim = dim(prevalence),
    dimnames = dimnames(prevalence)
  )

  list(prevalence = prevalence, deltas = deltas, incidence = incidence, pars = pars)

}

## define ODEs for SEIR
age_SEIR <- function(time, state, pars) {

  ## Each element is number of persons/total population (sum(state)=1)
  nr = nrow(pars$beta_I_c) ## Gettin the number of age groups, N
  state = matrix(
    state, nrow=nr,
    dimnames=list(
      paste0("AGE", 1:nr),
      c("S_u", "E_u", "C_u", "H_u", "R_u", "D_u", "S_v", "E_v", "C_v", "H_v", "R_v", "D_v")
    )
  )
  dstate = state * 0

  ## update rates in response to changes in hospitalisation
  pars %<>% add_rates(state)

  ## New E_u come from S being infected by E_u, C_u and H_u, E_v, C_v, H_v
  ## infectiousness reduction of vaccinated incorporated here
  new_E_u <- state[, "S_u"] * (
    state[, "E_u"] %*% pars$beta_E  +
    state[, "C_u"] %*% pars$beta_I_c +
    state[, "H_u"] %*% pars$beta_I_h +
    (1-pars$vax_infectious) * state[, "E_v"] %*% pars$beta_E +
    (1-pars$vax_infectious) * state[, "C_v"] %*% pars$beta_I_c +
    (1-pars$vax_infectious) * state[, "H_v"] %*% pars$beta_I_h
  )

  ## New E_v come from V being infected by E_u, C_u and H_u, E_v, C_v, H_v
  ## these are protected from infection by vax_infection
  ## infectiousness reduction of vaccinated incorporated here
  new_E_v <- state[, "S_v"] * (1-pars$vax_infection) * (
    state[, "E_u"] %*% pars$beta_E +
    state[, "C_u"] %*% pars$beta_I_c +
    state[, "H_u"] %*% pars$beta_I_h +
    (1-pars$vax_infectious) * state[, "E_v"] %*% pars$beta_E +
    (1-pars$vax_infectious) * state[, "C_v"] %*% pars$beta_I_c +
    (1-pars$vax_infectious) * state[, "H_v"] %*% pars$beta_I_h
  )

  ## Define absolute number of vaccinations not proportions (behaviour is not
  ## geometric) (i.e. not a proportion of S are vaccinated every day)
  vax <- state[, "S_u"]
  vax[] <- 0
  vaxleft <- pars$vax_rate

  ## Prioritised vaccination
  if(pars$vax_prioritised) {
    for(i in rev(seq_along(vax))) {
      if(vaxleft > state[i, "S_u"]) {
        vax[i] <- state[i, "S_u"]
        vaxleft <- vaxleft - state[i, "S_u"]
      } else {
        vax[i] <- vaxleft
        vaxleft <- 0
      }
    }
  } else {
    ## vaccinate a given proportion of all age groups equally
    if(pars$vax_rate > sum(state[, "S_u"])) vax <- state[, "S_u"]
    else vax <- pars$vax_rate*state[, "S_u"]/sum(state[, "S_u"])
  }

  ## S are lost to E_u by infection and to V by vaccination
  dstate[, "S_u"] = -new_E_u - vax

  ## E_u enter from S compartment, E_u leave as they enter sympomatic state in
  ## community or in hospital
  dstate[, "E_u"] = new_E_u - state[, "E_u"]*(pars$lambda_hu + pars$lambda_cu)

  ## C_u enter from E_u comparment, C_u leave as they die or recover in community
  dstate[, "C_u"] =
    state[, "E_u"]*pars$lambda_cu -
    state[, "C_u"]*(pars$sigma_cu + pars$mu_cu)

  ## H_u enter from E_u comparment, H_u leave as they die or recover in hospital
  dstate[, "H_u"] =
    state[, "E_u"]*pars$lambda_hu -
    state[, "H_u"]*(pars$sigma_hu + pars$mu_hu)

  ## R enter from community or hospital at respective recovery rates
  dstate[, "R_u"] =
    state[, "C_u"]*pars$sigma_cu +
    state[, "H_u"]*pars$sigma_hu

  ## D enter from community or hospital at respective mortality rates
  dstate[, "D_u"] =
    state[, "C_u"]*pars$mu_cu +
    state[, "H_u"]*pars$mu_hu

  ## S_v come from S_u by vaccination and are lost to E_v by infection
  dstate[, "S_v"] = vax - new_E_v

  ## E_v enter from V compartment, E_v leave as they enter sympomatic state in
  ## community or in hospital
  dstate[, "E_v"] = new_E_v - state[, "E_v"]*(pars$lambda_hv + pars$lambda_cv)

  ## I_v enter from E_v comparment, I_v leave as they die or recover in community
  dstate[, "C_v"] =
    state[, "E_v"]*pars$lambda_cv -
    state[, "C_v"]*(pars$sigma_cv + pars$mu_cv)

  ## I_v enter from E_v comparment, I_v leave as they die or recover in hospital
  dstate[, "H_v"] =
    state[, "E_v"]*pars$lambda_hv -
    state[, "H_v"]*(pars$sigma_hv + pars$mu_hv)

  ## R enter from community or hospital at respective recovery rates
  dstate[, "R_v"] =
    state[, "C_v"]*pars$sigma_cv +
    state[, "H_v"]*pars$sigma_hv

  ## D enter from community or hospital at respective mortality rates
  dstate[, "D_v"] =
    state[, "C_v"]*pars$mu_cv +
    state[, "H_v"]*pars$mu_hv

  list(dstate)

}

## calculate rates from periods
add_rates <- function(pars, state = NULL) {

  if(is.null(state)) hosp_required <- rep(0, length(pars$age_frac))
  else hosp_required <- state[, "H_u"] + state[, "H_v"]

  if(sum(hosp_required) > pars$hosp_capacity) {

    ## Define number that are actually hospitalised
    hospitalised <- hosp_required
    hospitalised[] <- 0
    hosp_left <- pars$hosp_capacity
    ## Fill in hospitalisations from oldest age group
    if(pars$hosp_prioritised) {
      for(i in rev(seq_along(hospitalised))) {
        if(hosp_left > hosp_required[i]) {
          hospitalised[i] <- hosp_required[i]
          hosp_left <- hosp_left - hospitalised[i]
        } else {
          hospitalised[i] <- hosp_left
          hosp_left <- 0
        }
      }
    } else {
      ## vaccinate a given proportion of all age groups equally
      hospitalised <- pars$hosp_capacity*hosp_required/sum(hosp_required)
    }

    ## total mortality of those that belong in hospital is weighted mean of
    ## proportion hospitalised and unhospitalised
    pars$total_hosp_mortality <- pars$hosp_mortality*hospitalised/hosp_required +
      pars$unhosp_mortality*(hosp_required - hospitalised)/hosp_required
    ## avoid dividing by zero when no hospitals required
    pars$total_hosp_mortality[hosp_required == 0] <- pars$hosp_mortality[hosp_required == 0]

  } else pars$total_hosp_mortality <- pars$hosp_mortality

  pars <- within(pars, {

    ##  E -> C_u (infectious, community, unvax)
    lambda_cu <- (1-prop_hosp)/incubation_period
    ##  E -> H_u (infectious, hospital, unvax)
    lambda_hu <- prop_hosp/incubation_period
    ##  E -> H_u (infectious, community, vax)
    lambda_cv <- (1-(1-vax_hosp)*prop_hosp)/incubation_period
    ##  E -> H_u (infectious, hospital, vax)
    lambda_hv <- (1-vax_hosp)*prop_hosp/incubation_period

    ## define rates from C_u -> R and D
    sigma_cu <- (1-comm_mortality)/comm_duration
    mu_cu <- comm_mortality/comm_duration

    ## define rates from C_v -> R and D
    sigma_cv <- (1-(1-vax_death)*comm_mortality)/comm_duration
    mu_cv <- (1-vax_death)*comm_mortality/comm_duration

    ## define rates from H_u -> R and D
    sigma_hu <- (1-total_hosp_mortality)/hosp_duration
    mu_hu <- total_hosp_mortality/hosp_duration

    ## define rates from H_v -> R and D
    sigma_hv <- (1-(1-vax_death)*total_hosp_mortality)/hosp_duration
    mu_hv <- (1-vax_death)*total_hosp_mortality/hosp_duration

  })

  return(pars)

}

## Converting from scenarios to betas using population structure
add_betas <- function(pars, mult=FALSE) {

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
        pars$scenario$isol_dur/pars$comm_duration *
        (1 - pars$scenario$isol_eff), 3)
  )

  ## proportion of contact rate maintained during stay in hospital
  ## home/school/work are removed, other are reduced by isol_eff
  hosp_coeffs <- c(0, 0, 0, 1 - pars$scenario$isol_eff)

  ## do social distancing and isolation have a combined effect?
  if (mult) isol_coeffs <- non_isol_coeffs * isol_coeffs
  else isol_coeffs <- pmin(non_isol_coeffs, isol_coeffs)

  ## calculation fraction of symptomatic cases that are isolated
  frac_isol <- pars$frac_symp *
    pars$scenario$isol *
    min(c(1, pars$scenario$isol_dur/pars$comm_duration))

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

  ## get total contact rates across different types
  polyscale_hosp <- map2(hosp_coeffs, pars$poly$scale[places], \(x, y) x*y) %>%
    Reduce("+", .)

  ## cocooning is a contact reduction specific to >65 year olds
  cocoon_coeffs = c(rep(1, 13), rep(1 - pars$scenario$cocoon_eff, 3))

  ## beta for Exposed cases takes:
  ## - contact rates unaffected by isolated but affected by social distancing
  ## - estimated transmission probability per contact
  ## - relative infectiousness of Exposed relative to Infected
  ## - cocooning effect on contact rates for elderly
  pars$beta_E = t(cocoon_coeffs * t(pars$inf_E_over_I * pars$p_trans * polyscale_non_isol))

  ## beta for Infected cases in the community takes:
  ## - weighted mean of contact rates from isolated and non-isolated Infecteds
  ## - estimated transmission probability per contact
  ## - cocooning effect on contact rates for elderly
  pars$beta_I_c = t(
    cocoon_coeffs *
    t(pars$p_trans *
      (frac_isol * polyscale_isol +
       (1 - frac_isol) * polyscale_non_isol))
  )

  ## beta for Infected cases in the hospital takes:
  ## - estimated transmission probability per contact
  ## - polyscale that removes home/work/school
  pars$beta_I_h = pars$p_trans * polyscale_hosp

  return(pars)

}

## get population data
## https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2022_Population1JanuaryByAge5GroupSex_Medium.zip
get_pop <- function(file = here("data/wpp_population.csv")) {
  fread(file)[Time %in% 2000:2024,
              .(ISO3_code, Location, Time, AgeGrp,
                AgeGrpStart, AgeGrpSpan, PopMale, PopFemale, PopTotal)] %>%
    as_tibble()
}


## get population fractions
get_pop_frac <- function(pop, iso3 = "DEU", yr = 2024) {

  pop %>%
    filter(ISO3_code == iso3 & Time == yr) %>%
    group_by(age_floor = ifelse(AgeGrpStart>=75, 75, AgeGrpStart)) %>%
    summarise(pop=1000*sum(PopTotal)) %>%
    ungroup() %>%
    mutate(pop_frac = prop.table(pop))

}

## compare different scenarios
vis_comparison <- function(scenarios,
                           what = c("prevalence", "deltas", "incidence"),
                           log = FALSE,
                           freescales = FALSE) {

  what <- match.arg(what)

  imap_dfr(scenarios, ~ mutate(extract(.x, what), scenario = .y)) %>%
    mutate(scenario = fct_inorder(scenario)) %>%
    group_by(day, compartment, scenario) %>%
    summarise(value = sum(value)) %>%
    ungroup() %>%
    ggplot(aes(day, value, color = scenario)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ compartment, scales = ifelse(freescales, "free_y", "fixed")) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0.01, 0.05)),
      trans = ifelse(log, "log10", "identity"),
      labels = scales::percent
    ) +
    scale_color_brewer(name = "Scenario", palette = "Dark2") +
    labs(x = "Day", y = "Proportion", color = "Category") +
    theme_minimal() +
    theme(
      legend.position = 'bottom',
      plot.background = element_rect(fill = "white")
    )

}

## compare different scenarios
vis_deaths <- function(scenarios,
                       what = c("prevalence", "deltas", "incidence"),
                       log = FALSE,
                       freescales = FALSE) {

  scen %>%
    map_dfc(
      ~ apply(.x$prevalence[365,,c("D_v", "D_u")], 1, sum)/pars$age_frac
    ) %>%
    mutate(age = fct_inorder(get_age_cat())) %>%
    pivot_longer(-age) %>%
    mutate(name = fct_inorder(name)) %>%
    ggplot(aes(age, value, color = name, group = name)) +
    geom_line() +
    geom_point() +
    scale_y_continuous(
      labels = scales::percent,
      expand = expansion(mult = c(0.01, 0.05))
    ) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      x = "Age category",
      y = "Deaths as proportion of population",
      color = "Scenario"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.background = element_rect("white")
    )

}

## show categories over time aggregated across age groups
vis_timeline <- function(out,
                         what = c("prevalence", "deltas", "incidence"),
                         log = FALSE,
                         freescales = FALSE) {

  what <- match.arg(what)
  apply(out[[what]], c(1, 3), sum) %>%
    {tibble(day = as.numeric(str_remove(rownames(.), "day_")), as_tibble(.))} %>%
    pivot_longer(-day, names_to = "compartment") %>%
    separate(compartment, c("compartment", "vax")) %>%
    mutate(
      vax = grepl("v", vax),
      compartment = fct_inorder(compartment)
    ) %>%
    ggplot(aes(day, value, color = compartment, linetype = vax)) +
    geom_line(linewidth = 2) +
    facet_wrap(~ compartment, scales = ifelse(freescales, "free_y", "fixed")) +
    scale_y_continuous(
      expand = expansion(mult = c(0.01, 0.05)),
      trans=ifelse(log, "log10", "identity")
    ) +
    scale_color_brewer(palette = "Dark2") +
    scale_linetype(name = "Vaccinated") +
    labs(x = "Day", y = "Proportion", color = "Category") +
    theme_minimal() +
    theme(legend.position = 'bottom')

}

## calculate various useful properties from transmission probability
get_R0 <- function(p_trans, pars) max(Re(eigen(get_matrices(p_trans, pars)$ngm)$values))
get_r <- function(p_trans, pars) max(Re(eigen(get_matrices(p_trans, pars)$delta)$values))
get_agesums <- function(p_trans, pars) max(Re(eigen(get_matrices(p_trans, pars)$agesums)$values))
get_doubling_time <- function(p_trans, pars) log(2)/get_r(p_trans, pars)
get_generation_time <- function(p_trans, pars) {
  mat <- get_matrices(p_trans, pars)
  max(Re(eigen(mat$agesums)$values))/max(Re(eigen(mat$ngm)$values))
}

## calculate different matrices used for R0 and r calculations
## these assume a completely susceptible population!!
get_matrices <- function(p_trans, pars) {

  nr <- length(pars$age_frac)

  ## these are new infections (E_u to E_u, C_u to E_u, I_h to E_u) (filling by column)
  T <- cbind(
    ## E generate E with contact_rate*p_transmission*infectiousness_ratio
    rbind(t(p_trans * pars$inf_E_over_I * pars$poly$mod$all), diag(0, nr), diag(0, nr)),
    ## I_c generate E with contact_rate*p_transmission
    rbind(t(p_trans*pars$poly$mod$all), diag(0, nr), diag(0, nr)),
    ## I_h generate E with contact_rate*p_transmission using ONLY other contacts
    rbind(t(p_trans*pars$poly$mod$other), diag(0, nr), diag(0, nr))
  )

  pars$age_frac * t(pars$poly$scale$all)

  ## these are the transitions (filling by column) of individuals between
  ## categories (lots of zeroes because no transitions between age-groups)
  Eps <- cbind(
    rbind(-diag(pars$lambda_cu + pars$lambda_hu, nr),
          diag(pars$lambda_cu, nr),
          diag(pars$lambda_hu, nr)),
    rbind(diag(0, nr), -diag(pars$sigma_cu + pars$mu_cu, nr), diag(0, nr)),
    rbind(diag(0, nr), diag(0, nr), -diag(pars$sigma_hu + pars$mu_hu, nr))
  )

  ## average future time spent by index j in index i
  Eps_inv <- cbind(
    rbind(
      ## E will spend in E given by inverse of the rates leave into I_h and I_c
      diag(1/(pars$lambda_cu + pars$lambda_hu), nr),
      ## average future time spent by E in I_c: the first part is the duration
      ## I_c spends in the community, the second part is the proportion of E
      ## that go into I-c
      diag(
        1/(pars$sigma_cu + pars$mu_cu) * pars$lambda_cu/(pars$lambda_cu + pars$lambda_hu),
        nr
      ),
      ## average future time spent by E in I_h
      diag(
        1/(pars$sigma_hu + pars$mu_hu)*pars$lambda_hu/(pars$lambda_cu + pars$lambda_hu),
        nr
      )
    ),
    ## I_c spends comm_duration in community
    rbind(diag(0, nr), diag(1/(pars$sigma_cu + pars$mu_cu), nr), diag(0, nr)),
    ## I_h spends hosp_duration in community
    rbind(diag(0, nr), diag(0, nr), diag(1/(pars$sigma_hu + pars$mu_hu), nr))
  )

  ## calculate the mean "age" of infections caused by each compartment
  age_E <- if(length(pars$incubation_period) == 1) rep(pars$incubation_period, nr)/2
           else pars$incubation_period/2
  age_C <- pars$incubation_period + if(length(pars$comm_duration) == 1) rep(pars$comm_duration, nr)/2
                                    else pars$comm_duration/2
  age_H <- pars$incubation_period + if(length(pars$hosp_duration) == 1) rep(pars$hosp_duration, nr)/2
                                    else pars$hosp_duration/2
  ages <- c(age_E, age_C, age_H)

  ## The RHS multiplies the days spent by j in i times the average age of the
  ## case they would generate. The matrix multiplication operator then
  ## incorporates the force of infection (i.e. daily infections caused): if j
  ## spends 5 days in i, the average age caused is 7 days, and infects on
  ## average 2 cases per day, the total "age days" generated by j in i in 5*7*2
  ## = 70. The matrix multiplication operator then sums across the infector j to
  ## give the total "age days" it contributes to compartment i over the course
  ## of it's lifetime.
  agesums <- T %*% apply(Eps_inv, 2, \(x) x * ages)

  ## This looks just at change and is used for r
  delta <- T + Eps

  ## This gives new infections and is used for R0
  ngm <- T %*% Eps_inv

  list(T = T, Eps = Eps, Eps_inv = Eps_inv, delta = delta, ngm = ngm, agesums = agesums)

}

## convert transmission probability to doubling time in absence of any
## interventions
trans_prob_to_doubling_time <- function(p_trans, pars) {

  ## E infects E according to contact rate and probability and leaves its own
  ## compartment according to lambda rates (into community and hosp)
  nr <- length(pars$age_frac)
  E_to_E <-
    t(p_trans*pars$inf_E_over_I*pars$poly$mod$all) -
    diag(pars$lambda_c + pars$lambda_h, nr)
  E_to_I_c <- diag(pars$lambda_c, nr)
  E_to_I_h <- diag(pars$lambda_h, nr)

  I_c_to_E <- t(p_trans*pars$poly$mod$all)
  I_c_to_I_c <- -diag(pars$sigma_c + pars$mu_c, nr)
  I_c_to_I_h <- diag(0, nr)

  ## this is assuming that in the absence of interventions (for calculating R0),
  ## hospitalised cases only have contacts of type "other"
  I_h_to_E <- t(p_trans*pars$poly$mod$all)
  I_h_to_I_c <- diag(0, nr)
  I_h_to_I_h <- -diag(pars$sigma_h + pars$mu_h, nr)

  b <- cbind(
    rbind(E_to_E, E_to_I_c, E_to_I_h),
    rbind(I_c_to_E, I_c_to_I_c, I_c_to_I_h),
    rbind(I_h_to_E, I_h_to_I_c, I_h_to_I_h)
  )

  res <- log(2)/max(Re(eigen(b)$values))
  ifelse(res<0, Inf, res)

}

## Doubling time and principal eigenvector from next generation matrix
next_gen_analysis <- function(pars, verbose=FALSE) {
  with(pars, {
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

## list of potential scenarios and their impact on contact rates (in home,
## school, work other) as well as the proportion of cases isolated (isol), the
## effectiveness of this isolation (isol_eff) and the duration of isolaion
## (isol_dur)
get_scenarios <- function() {

readHTMLTable(here("notes/intervention_scenarios.html"))[[1]] %>%
  mutate(across(-Scenario, as.numeric)) %>%
  split(.$Scenario) %>%
  map(~ discard_at(as.list(.x), "Scenario"))

}

## extract data from solved ODE and convert to DF
extract <- function(solved,
                    what = c("prevalence", "deltas", "incidence")) {

  what <- match.arg(what)
  apply(solved[[what]], c(1, 3), sum) %>%
    {tibble(day = as.numeric(str_remove(rownames(.), "day_")), as_tibble(.))} %>%
    pivot_longer(-day, names_to = "compartment") %>%
    separate(compartment, c("compartment", "vax")) %>%
    mutate(
      vax = grepl("v", vax),
      compartment = fct_inorder(compartment)
    )

}

## ifr calculations from https://link.springer.com/article/10.1007/s10654-020-00698-1#Sec7
age_to_ifr <- function(age) 10^(-3.27 + 0.0524*age)/100

## get age categories
get_age_cat <- function() {
  x <- map_chr(age$age_floor, ~ paste0(.x, "-", .x+4))
  x[length(x)] <- "75+"
  return(x)
}

## plot saving function
save_plot <- function(p, file,
                      width = 20.16, height = 13.35,
                      units = 'cm',
                      dpi = 600,
                      folder = "figures",
                      ...) {

  if(is.null(p)) return(NULL)

  if(!file.exists(here(folder))) dir.create(here(folder))

  ggsave(
    filename = here(folder, file),
    plot = p,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    ...
    )

}
