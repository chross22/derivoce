#' Current speed from velocity components
#'
#' The magnitude of the horizontal velocity vector, `sqrt(u^2 + v^2)`. Reproduces
#' the older pipeline's `uv`. Pass the result through [horizontal_gradient()] to
#' get its spatial gradient, the older pipeline's `uv_grad`.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`
#' @param u name of the eastward velocity column
#' @param v name of the northward velocity column
#' @param name name for the new column
#' @return `env_dat` with a current speed column added, in the units of `u`/`v`
#' @examples
#' \dontrun{
#' env <- current_speed(env, "uo", "vo") |> horizontal_gradient("speed")
#' }
#' @export
current_speed <- function(env_dat, u = "UO", v = "VO", name = "speed") {
  resolve_vars(env_dat, c(u, v))
  env_dat[[name]] <- sqrt(env_dat[[u]]^2 + env_dat[[v]]^2)
  env_dat
}

#' Eddy kinetic energy
#'
#' The kinetic energy carried by the *departure* of the flow from its mean:
#'
#' \deqn{EKE = \tfrac{1}{2}(u'^2 + v'^2), \quad u' = u - \bar{u}}
#'
#' High EKE marks an energetic, variable region — meanders, rings, and eddies —
#' as distinct from a strong but steady current, which carries high total kinetic
#' energy but little eddy energy.
#'
#' @section What the anomaly is measured against:
#' EKE is only defined relative to a mean, and the choice of mean decides what
#' counts as an "eddy". This is a scientific decision, not an implementation
#' detail, so it is an argument with no silently-correct default:
#'
#' \itemize{
#'   \item `"record"` (the default) subtracts the mean over the whole time series
#'     at each location. This is the textbook definition. Because the seasonal
#'     cycle of the mean circulation is left in the anomaly, a region whose
#'     currents merely strengthen every summer will register as energetic.
#'   \item `"climatology"` subtracts a separate mean for each calendar month, so
#'     the repeatable seasonal cycle is removed and only departures *from the
#'     usual conditions for that month* count. Use this when the question is about
#'     anomalous years rather than about which places are energetic.
#'   \item An integer subtracts a centred rolling mean over that many time steps,
#'     so slow drift in the mean state is removed and only variability faster than
#'     the window survives.
#' }
#'
#' Means are computed per location, so a spatially varying mean circulation is
#' removed correctly rather than being smeared into the anomaly.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`
#' @param u name of the eastward velocity column, in m/s
#' @param v name of the northward velocity column, in m/s
#' @param reference `"record"`, `"climatology"`, or a positive integer number of
#'   time steps for a centred rolling mean
#' @param name name for the new column
#' @return `env_dat` with an EKE column added, in m^2/s^2
#' @examples
#' \dontrun{
#' env <- eke(env, "uo", "vo")                          # against the record mean
#' env <- eke(env, "uo", "vo", reference = "climatology") # seasonal cycle removed
#' env <- eke(env, "uo", "vo", reference = 12)          # 12-step rolling mean
#' }
#' @export
eke <- function(env_dat, u = "UO", v = "VO", reference = "record", name = "EKE") {
  resolve_vars(env_dat, c(u, v))

  if (is.numeric(reference) && reference < 2) {
    stop("A numeric reference must be at least 2; a rolling mean of one step ",
         "would make every anomaly zero.", call. = FALSE)
  }
  if (is.character(reference) && !(reference %in% c("record", "climatology"))) {
    stop("reference must be \"record\", \"climatology\", or a positive integer.",
         call. = FALSE)
  }

  u_anomaly <- anomaly(env_dat, u, reference)
  v_anomaly <- anomaly(env_dat, v, reference)

  env_dat[[name]] <- 0.5 * (u_anomaly^2 + v_anomaly^2)
  env_dat
}

#' Departure of a covariate from its mean, per location
#'
#' @param env_dat an `sf` POINT object
#' @param var covariate column name
#' @param reference `"record"`, `"climatology"`, or a rolling window length
#' @return numeric vector of anomalies, one per row
#' @keywords internal
anomaly <- function(env_dat, var, reference) {
  values <- env_dat[[var]]
  location <- location_key(env_dat)

  group <- switch(
    as.character(reference)[1],
    record = location,
    # One mean per location AND calendar month, so the seasonal cycle is removed.
    climatology = paste(location, env_dat$MONTH, sep = "|"),
    NULL
  )

  if (!is.null(group)) {
    means <- tapply(values, group, mean, na.rm = TRUE)
    return(values - as.numeric(means[group]))
  }

  rolling_anomaly(env_dat, values, location, window = reference)
}

#' Anomaly against a centred rolling mean
#'
#' @param env_dat an `sf` POINT object
#' @param values the covariate values
#' @param location per-row location key
#' @param window number of time steps in the rolling mean
#' @return numeric vector of anomalies
#' @keywords internal
rolling_anomaly <- function(env_dat, values, location, window) {
  steps <- time_steps(env_dat)
  step_index <- match_step_index(env_dat, steps)
  half <- (window - 1) / 2

  result <- rep(NA_real_, length(values))

  for (i in seq_len(nrow(steps))) {
    current <- which(step_index == i)
    # Centred, and clipped at the ends of the record rather than returning NA
    # there - a shorter window at the edges is more useful than no value.
    contributing <- seq(max(1, i - ceiling(half)),
                        min(nrow(steps), i + floor(half)))

    total <- rep(0, length(current))
    counted <- rep(0, length(current))
    for (j in contributing) {
      other <- which(step_index == j)
      matched <- values[other][match(location[current], location[other])]
      total <- total + ifelse(is.na(matched), 0, matched)
      counted <- counted + !is.na(matched)
    }
    means <- ifelse(counted > 0, total / counted, NA_real_)
    result[current] <- values[current] - means
  }
  result
}
