#' Lagged covariate values
#'
#' Adds each covariate's value from `n` units earlier at the same location.
#' Populations respond to conditions with a delay. A bloom feeds the animals
#' sampled a month later, not the ones sampled during it, so the lagged value can
#' carry more signal than the concurrent one.
#'
#' @section Lag by calendar units:
#' `by` is required and says what `n` counts: `"hour"`, `"day"`, `"month"` or
#' `"year"`. `n = 3` with `by = "month"` finds the step stamped exactly three calendar
#' months earlier, and returns `NA` if there is no such step.
#'
#' There is no default and no option to count positions in the series. "Three
#' months ago" is a statement about the organism; "three steps ago" is a
#' statement about how the data happened to be fetched, and the two stop agreeing
#' the moment a month is missing. In a monthly series missing April, a position
#' count would treat March as May's predecessor and quietly turn a one-month lag
#' into a two-month one. A calendar lag returns `NA` for May instead, because
#' April genuinely is not there. With no default, the unit is chosen on purpose:
#' a guess would be right for a daily record and wrong for a monthly one.
#'
#' Calendar lags are matched on the exact `YEAR`/`MONTH`/`DAY` stamp. For monthly
#' products, whose day is always 1, that is exact. For daily products,
#' `by = "month"` from the 31st looks for a 31st, and months that have no 31st
#' return `NA` rather than silently sliding to the 30th.
#'
#' Locations are matched by coordinate, so this assumes a fixed grid across time
#' steps, which is what gridded products give.
#'
#' @section Reproducing the published lag:
#' Ross et al. (2023) used a **one-month** lag of sea surface temperature:
#'
#' ```
#' lag_covariate(env, "SST", n = 1, by = "month")
#' ```
#'
#' @section Several lags at once:
#' `n` may be a vector, which adds one column per lag. This is what an
#' autoregressive design needs: a covariate's own past at a series of offsets,
#' entered together as predictors.
#'
#' ```
#' lag_covariate(env, "SST", n = 1:3, by = "year")
#' # adds SST_lag1year, SST_lag2year, SST_lag3year
#' ```
#'
#' With `by = "year"` this gives the same calendar month in each preceding year,
#' so the seasonal cycle is held fixed and what remains is the interannual
#' signal.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param vars covariate columns to lag; `NULL` does all of them
#' @param n how far to look back, counted in `by` units. A vector adds one
#'   column per lag.
#' @param by what `n` counts: `"hour"`, `"day"`, `"month"` or `"year"`.
#'   Required.
#' @param suffix suffix for the new columns. The default names both `n` and the
#'   unit, as in `SST_lag1month`. Supply one entry per lag when `n` has several.
#' @return `env_dat` with a lagged column per covariate and lag. Steps with no
#'   predecessor are `NA`.
#' @references
#' Ross C, Runge J, Roberts J, Brady D, Tupper B, Record N (2023). Estimating
#' North Atlantic right whale prey based on Calanus finmarchicus thresholds.
#' *Marine Ecology Progress Series* **703**, 1-16. \doi{10.3354/meps14204}
#' @examples
#' \dontrun{
#' env <- lag_covariate(env, "SST", by = "month")        # SST_lag1month
#' env <- lag_covariate(env, "CHL", n = 3, by = "month") # CHL_lag3month
#' env <- lag_covariate(env, "SST", n = 1, by = "year")  # same month last year
#' env <- lag_covariate(env, "SST", n = 30, by = "day")  # daily products
#'
#' # An autoregressive set: the same month in each of the last three years
#' env <- lag_covariate(env, "SST", n = 1:3, by = "year")
#' }
#' @export
lag_covariate <- function(env_dat, vars = NULL, n = 1, by, suffix = NULL) {
  by <- calendar_unit(if (missing(by)) NULL else by, "lag_covariate")
  vars <- resolve_vars(env_dat, vars, kind = "temporal")
  if (!is.numeric(n) || length(n) < 1 || anyNA(n) || any(n < 1) ||
      any(n != round(n))) {
    stop("n must be one or more whole numbers of at least 1.", call. = FALSE)
  }
  if (!is.null(suffix) && length(n) > 1 && length(suffix) != length(n)) {
    stop("`suffix` must have one entry per lag when `n` has several, or be ",
         "NULL to name them automatically.", call. = FALSE)
  }
  suffixes <- suffix %||% paste0("_lag", n, by)

  steps <- time_steps(env_dat)
  location <- location_key(env_dat)
  step_index <- match_step_index(env_dat, steps)

  for (k in seq_along(n)) {
    source_step <- lag_source_step(steps, n[k], by)

    for (v in vars) {
      lagged <- rep(NA_real_, nrow(env_dat))

      for (i in seq_len(nrow(steps))) {
        from <- source_step[i]
        if (is.na(from)) next
        current <- which(step_index == i)
        earlier <- which(step_index == from)

        # Match by location rather than assuming the two steps list their
        # points in the same order.
        lagged[current] <- env_dat[[v]][earlier][match(location[current],
                                                       location[earlier])]
      }
      env_dat[[paste0(v, suffixes[k])]] <- lagged
    }
  }
  env_dat
}

#' The calendar unit a time function counts in
#'
#' Every function that looks back in time counts in calendar units, never in
#' positions in the series, and none has a default unit. Positions are what a gap
#' silently corrupts, and a default unit would be right for one cadence and wrong
#' for another, so the caller must say which.
#'
#' @param by the value the caller supplied, or `NULL` if none was given
#' @param fn name of the calling function, for the message
#' @return `"hour"`, `"day"`, `"month"` or `"year"`
#' @keywords internal
calendar_unit <- function(by, fn) {
  units <- c("hour", "day", "month", "year")
  if (is.null(by)) {
    stop(fn, "() needs `by`, the calendar unit to count in: \"hour\", ",
         "\"day\", \"month\" or \"year\".",
         "\n  There is no default, because the right unit depends on the data ",
         "(a daily record and a monthly one differ) and a wrong guess would be ",
         "silent.", call. = FALSE)
  }
  if (identical(by, "step")) {
    stop(fn, "() no longer counts positions in the series (`by = \"step\"`).",
         "\n  Across a gap a position count reaches the wrong month without ",
         "saying so. Use a calendar unit: \"hour\", \"day\", \"month\" or ",
         "\"year\".",
         call. = FALSE)
  }
  if (!is.character(by) || length(by) != 1 || !by %in% units) {
    stop(fn, "(): `by` must be one of \"hour\", \"day\", \"month\" or \"year\".",
         call. = FALSE)
  }
  by
}

#' Which earlier step each step draws its lagged value from
#'
#' Returns one index per time step, or `NA` where the lag lands outside the
#' record. Separating this from [lag_covariate()] keeps the calendar arithmetic
#' in one place and lets it be tested on the step table alone.
#'
#' @param steps a time-step table from [time_steps()]
#' @param n how far back, in `by` units
#' @param by `"hour"`, `"day"`, `"month"`, or `"year"`
#' @return integer vector of source step indices, `NA` where there is none
#' @keywords internal
lag_source_step <- function(steps, n, by) {
  # On a sub-daily record a calendar lag means the same hour, one day or one
  # month before; without the hour in the key, every hour of a day would draw
  # from whichever hour happened to be listed first.
  hour <- if ("HOUR" %in% names(steps)) steps$HOUR else 0L

  if (identical(by, "hour")) {
    # Whole hours since a fixed origin, so a lag crosses midnight, month ends and
    # leap days without any calendar reasoning.
    current <- as.numeric(as.Date(paste(steps$YEAR, steps$MONTH, steps$DAY,
                                        sep = "-"))) * 24 + hour
    return(match(current - n, current))
  }

  if (identical(by, "day")) {
    # Real dates, so month lengths and leap years take care of themselves.
    current <- as.Date(paste(steps$YEAR, steps$MONTH, steps$DAY, sep = "-"))
    return(match(paste(current - n, hour), paste(current, hour)))
  }

  # Months and years are the same arithmetic, done on a month counter so that
  # subtracting never has to reason about how long a month is.
  months_back <- if (identical(by, "year")) n * 12 else n
  counter <- steps$YEAR * 12 + (steps$MONTH - 1)
  target <- counter - months_back

  match(paste(target %/% 12, target %% 12 + 1, steps$DAY, hour),
        paste(steps$YEAR, steps$MONTH, steps$DAY, hour))
}

#' Time-integrated covariate
#'
#' Accumulates a covariate over preceding time steps at each location. A survey
#' does not sample the food available that instant so much as the food that has
#' built up since the season began, which is what the original pipeline's
#' `int_chl` captured for chlorophyll.
#'
#' The default `window = "year"` reproduces `int_chl` of Ross et al. (2023),
#' which integrated chlorophyll from January: a running sum from January
#' of each year, reset at the year boundary. A numeric window instead sums over
#' that many trailing time steps, giving a rolling total that does not reset.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param vars covariate columns to integrate; `NULL` does all of them
#' @param window `"year"` to accumulate from the start of each calendar year,
#'   `"all"` to accumulate over the whole record, or a positive integer for a
#'   rolling window of that many time steps
#' @param suffix suffix for the new columns
#' @return `env_dat` with an integrated column per covariate
#' @references
#' Ross C, Runge J, Roberts J, Brady D, Tupper B, Record N (2023). Estimating
#' North Atlantic right whale prey based on Calanus finmarchicus thresholds.
#' *Marine Ecology Progress Series* **703**, 1-16. \doi{10.3354/meps14204}
#' @examples
#' \dontrun{
#' env <- integrate_covariate(env, "CHL")               # int_chl, Ross et al. 2023
#' env <- integrate_covariate(env, "CHL", window = 3)   # trailing 3-step total
#' }
#' @export
integrate_covariate <- function(env_dat, vars = NULL, window = "year",
                                suffix = "_int") {
  vars <- resolve_vars(env_dat, vars, kind = "temporal")
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
