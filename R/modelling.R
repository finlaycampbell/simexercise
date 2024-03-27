## load packages
pacman::p_load(XML, here, data.table, countrycode, readxl, glue, tidyverse, rio, deSolve, magrittr)

## socialmixr?
## quarto dashboard?

## source functions
source(here("R/functions.R"))

## get various scenarios defining intervention effects
scenarios <- get_scenarios()

## load country data
cdat <- import(here("data/cdat.rda"))

## get countrynames matched to ISO3 code for polymod dataset
countrynames <- get_countrynames()

pars <- get_pars()
solved <- solve_ode(pars)
vis_timeline(solved, "prevalence", log = FALSE, freescales = TRUE)

test_vaccine() %>% save_plot("vaccine_deaths.png")

test_hosp() %>% save_plot("hosp_deaths.png")

scen <- list(
  "No vaccination" = list(vax_rate = 0, vax_prioritised = FALSE),
  "Prioritised vaccination for elderly" = list(vax_rate = 0.002, vax_prioritised = TRUE),
  "Equal vaccination for all ages" = list(vax_rate = 0.002, vax_prioritised = FALSE)
) %>%
  map(~ solve_ode(get_pars(vax_rate = .x[[1]], vax_prioritised = .x[[2]])))

scen %>%
  vis_comparison(freescales = TRUE) %>%
  save_plot("vaccine_trajectories.png")

test_vaccine <- function() {

  lbs <- c(
    novax = "No vaccination",
    vaxunpriori = "Equal vaccination for all ages",
    vaxpriori = "Prioritised vaccination for elderly"
  )

  df <- list(
    novax = list(vax_rate = 0, vax_prioritised = FALSE),
    vaxpriori = list(vax_rate = 0.002, vax_prioritised = TRUE),
    vaxunpriori = list(vax_rate = 0.002, vax_prioritised = FALSE)
  ) %>%
    map(~ runODE(get_pars(R0 = 5, incubation_period = 3, generation_time = 9,
                          vax_rate = .x[[1]], vax_prioritised = .x[[2]]))) %>%
    map_dfc(
      ~ apply(.x$prevalence[365,,c("D_v", "D_u")], 1, sum)/pars$age_frac
    )

  df %>%
    mutate(age = fct_inorder(get_age_cat())) %>%
    pivot_longer(-age) %>%
    mutate(name = factor(name, names(lbs))) %>%
    ggplot(aes(age, value, color = name, group = name)) +
    geom_line() +
    geom_point() +
    scale_y_continuous(
      labels = scales::percent,
      expand = expansion(mult = c(0.01, 0.05))
    ) +
    scale_color_brewer(palette = "Dark2", labels = lbs) +
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

test_hosp <- function() {

  df <- list(
    hosp1 = list(hosp_capacity = 0.000, hosp_prioritised = FALSE),
    ## hosp2 = list(hosp_capacity = 0.001, hosp_prioritised = FALSE),
    hosp3 = list(hosp_capacity = 0.001, hosp_prioritised = TRUE),
    hosp3 = list(hosp_capacity = 0.002, hosp_prioritised = TRUE),
    hosp4 = list(hosp_capacity = 0.003, hosp_prioritised = TRUE),
    hosp4 = list(hosp_capacity = 0.005, hosp_prioritised = TRUE),
    hosp5 = list(hosp_capacity = 0.008, hosp_prioritised = TRUE)
  ) %>%
    map(~ runODE(get_pars(vax_rate = 0, hosp_capacity = .x[[1]], hosp_prioritised = .x[[2]]))) %>%
    map_dfc(
      ~ apply(.x$prevalence[365,,c("D_v", "D_u")], 1, sum)#/pars$age_frac
    )

  df %>%
    mutate(age = seq_len(n())) %>%
    pivot_longer(-age) %>%
    ggplot(aes(age, value, color = name, group = name)) +
    geom_line() +
    geom_point() +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal()

}
