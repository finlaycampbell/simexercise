## get polymod contact matrices (POLYMOD + Prem et al.))
## https://doi.org/10.1371/journal.pcbi.1005697.s002
import_polymod <- function(iso3 = "DEU", yr = 2024,
                           loc = c("all_locations", "home", "other_locations",
                                   "school", "work")) {

  loc <- match.arg(loc)
  as.matrix(
    read_excel(
      paste0(
        system.file("extdata/polymod", "my_raw_data.csv", package="simex"),
        loc, "_", countrynames[COUNTRY_CODE==iso3, Set], ".xlsx"),
      sheet=countrynames[COUNTRY_CODE==iso3, COUNTRY_NAME],
      col_names=(countrynames[COUNTRY_CODE==iso3, Set]==1),
      .name_repair=toupper
    )
  )

}

## get all polymod data in one list
get_poly <- function(iso3 = "DEU", yr = 2024,
                     file = here("data/polymod/MUestimates_")) {

  poly <- list()
  poly$mod <- map(
    c("all_locations", "home", "school", "work", "other_locations"),
    import_polymod, iso3 = iso3, yr = yr, file = file
  ) %>%
    set_names(c("all", "home", "school", "work", "other"))

  ## get population sizes
  populations <- get_pop_frac(pop, iso3, yr = 2024)$pop

  ## Regularize the darned polymod matrix for the national population distribution
  ## Arithmetic mean, to ensure additivity and consistency with flu-models approach
  poly$mod %<>% map(~ ((populations * .x) + t(populations * .x))/(2 * populations))

  ## Scale the polymods
  poly$scale <- map(poly$mod, ~ t(t(.x)/prop.table(populations)))

  return(poly)

}
