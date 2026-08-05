#' Lagged covariate values
#'
#' Adds each covariate's value `n` time steps earlier at the same location.
#' Populations respond to conditions with a delay — a bloom feeds the animals
#' sampled a month later, not the ones sampled during it — so the lagged value can
#' carry more signal than the concurrent one.
#'
#' The default `n = 1` reproduces the original pipeline's `lag_sst`.
#'
#' Locations are matched by coordinate, so this assumes a fixed grid across time
#' steps, which is what gridded products give.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`
#' @param vars covariate columns to lag; `NULL` does all of them
#' @param n number of time steps to look back
#' @param suffix suffix for the new columns; the default includes `n`
#' @return `env_dat` with a lagged column per covariate. The first `n` time steps
#'   have no predecessor and are `NA`.
#' @examples
#' \dontrun{
#' env <- lag_covariate(env, "SST")           # SST_lag1, the original's lag_sst
#' env <- lag_covariate(env, "CHL", n = 2)    # CHL_lag2
#' }
#' @export
lag_covariate <- function(env_dat, vars = NULL, n = 1, suffix = NULL) {
  vars <- resolve_vars(env_dat, vars)
  if (n < 1) stop("n must be at least 1.", call. = FALSE)
  suffix <- suffix %||% paste0("_lag", n)

  steps <- time_steps(env_dat)
  location <- location_key(env_dat)
  step_index <- match_step_index(env_dat, steps)

  for (v in vars) {
    lagged <- rep(NA_real_, nrow(env_dat))

    for (i in seq_len(nrow(steps))) {
      if (i - n < 1) next
      current <- which(step_index == i)
      earlier <- which(step_index == i - n)

      # Match by location rather than assuming the two steps list their points
      # in the same order.
      lagged[current] <- env_dat[[v]][earlier][match(location[current],
                                                     location[earlier])]
    }
    env_dat[[paste0(v, suffix)]] <- lagged
  }
  env_dat
}

#' Time-integrated covariate
#'
#' Accumulates a covariate over preceding time steps at each location. A survey
#' does not sample the food available that instant so much as the food that has
#' built up since the season began, which is what the original pipeline's
#' `int_chl` captured for chlorophyll.
#'
#' The default `window = "year"` reproduces `int_chl`: a running sum from January
#' of each year, reset at the year boundary. A numeric window instead sums over
#' that many trailing time steps, giving a rolling total that does not reset.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`
#' @param vars covariate columns to integrate; `NULL` does all of them
#' @param window `"year"` to accumulate from the start of each calendar year,
#'   `"all"` to accumulate over the whole record, or a positive integer for a
#'   rolling window of that many time steps
#' @param suffix suffix for the new columns
#' @return `env_dat` with an integrated column per covariate
#' @examples
#' \dontrun{
#' env <- integrate_covariate(env, "CHL")               # the original's int_chl
#' env <- integrate_covariate(env, "CHL", window = 3)   # trailing 3-step total
#' }
#' @export
integrate_covariate <- function(env_dat, vars = NULL, window = "year",
                                suffix = "_int") {
  vars <- resolve_vars(env_dat, vars)
  if (is.numeric(window) && window < 1) {
    stop("A numeric window must be at least 1.", call. = FALSE)
  }
  if (is.character(window) && !(window %in% c("year", "all"))) {
    stop("window must be \"year\", \"all\", or a positive integer.", call. = FALSE)
  }

  steps <- time_steps(env_dat)
  location <- location_key(env_dat)
  step_index <- match_step_index(env_dat, steps)

  for (v in vars) {
    integrated <- rep(NA_real_, nrow(env_dat))

    for (i in seq_len(nrow(steps))) {
      contributing <- contributing_steps(steps, i, window)
      current <- which(step_index == i)

      total <- rep(0, length(current))
      counted <- rep(0, length(current))
      for (j in contributing) {
        earlier <- which(step_index == j)
        values <- env_dat[[v]][earlier][match(location[current], location[earlier])]
        # NA means that location is absent from that step, not a zero
        # contribution, so track how many steps actually contributed.
        total <- total + ifelse(is.na(values), 0, values)
        counted <- counted + !is.na(values)
      }
      total[counted == 0] <- NA_real_
      integrated[current] <- total
    }
    env_dat[[paste0(v, suffix)]] <- integrated
  }
  env_dat
}

#' Time steps contributing to an integration window
#'
#' @param steps a time-step table from `time_steps()`
#' @param i index of the step being computed
#' @param window `"year"`, `"all"`, or a number of trailing steps
#' @return integer indices into `steps`, inclusive of `i`
#' @keywords internal
contributing_steps <- function(steps, i, window) {
  if (identical(window, "all")) return(seq_len(i))
  if (identical(window, "year")) {
    return(which(steps$YEAR == steps$YEAR[i])[which(steps$YEAR == steps$YEAR[i]) <= i])
  }
  seq(max(1, i - window + 1), i)
}

#' A stable key for each point's location
#'
#' Rounded so that coordinates written and re-read at different precisions still
#' match across time steps.
#'
#' @param env_dat an `sf` POINT object
#' @param digits rounding applied before forming the key
#' @return character vector, one per row
#' @keywords internal
location_key <- function(env_dat, digits = 6) {
  coords <- round(sf::st_coordinates(env_dat), digits)
  paste(coords[, 1], coords[, 2], sep = ",")
}
