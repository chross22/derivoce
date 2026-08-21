#' Potential density of seawater
#'
#' Density is what stratification, mixing and buoyancy actually depend on, and
#' temperature alone is a poor stand-in for it wherever salinity varies. The
#' Gulf of Maine is one of those places: Scotian Shelf inflow arrives cold *and*
#' fresh, and the two pull the density in opposite directions, so a cold anomaly
#' can be either denser or lighter than the water it displaces depending on how
#' fresh it is.
#'
#' Uses the UNESCO (1983) equation of state at one atmosphere. Copernicus
#' `thetao` is already a potential temperature, so applying it here gives
#' potential density directly, conventionally reported as sigma-theta: density
#' in kg/m^3 minus 1000.
#'
#' "One atmosphere" describes the *pressure* the density is referenced to, not
#' the depth the temperature and salinity came from. Potential density is
#' exactly the quantity you want for water sampled at depth — it is what that
#' water would weigh if brought to the surface, which is what makes two levels
#' comparable. `datamatch::accessCopernicus()` takes a `depth` argument, so this
#' applies to any level, and [buoyancy_frequency()] uses two of them.
#'
#' @section What it is not:
#' This is the density a parcel would have if brought to the surface. It is the
#' right quantity for comparing water masses and for deciding what floats over
#' what, and it deliberately ignores pressure, so it is **not** in-situ density
#' and should not be used where the compressibility of deep water matters.
#'
#' Applying it to a temperature that is not a potential temperature gives
#' in-situ surface density instead, which is the same number at the surface and
#' increasingly wrong with depth.
#'
#' @section Range:
#' The polynomial is fitted over roughly -2 to 40 degrees C and 0 to 42 PSU.
#' Outside that it still returns a number, and that number is an extrapolation
#' of a fit rather than a density, so values beyond the range are warned about.
#' Fresh water from a river mouth and ice-melt surface layers are the usual
#' culprits.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param temperature name of the temperature column, in degrees C. Should be a
#'   potential temperature
#' @param salinity name of the salinity column, in PSU
#' @param sigma return sigma-theta, that is density minus 1000. `FALSE` returns
#'   absolute density
#' @param name name for the new column
#' @return `env_dat` with a density column, in kg/m^3
#' @references
#' UNESCO (1983). Algorithms for computation of fundamental properties of
#' seawater. *UNESCO Technical Papers in Marine Science* **44**.
#' @examples
#' \dontrun{
#' env <- potential_density(env)
#'
#' # Density anomalies say more about inflow than temperature anomalies do,
#' # because they combine the cold and the fresh into one number.
#' env <- cell_anomaly(env, "sigma_theta")
#' }
#' @seealso [water_mass_fraction()], [vertical_gradient()], [cell_anomaly()]
#' @export
potential_density <- function(env_dat, temperature = "SST", salinity = "SSS",
                              sigma = TRUE, name = NULL) {
  resolve_vars(env_dat, c(temperature, salinity))

  t <- env_dat[[temperature]]
  s <- env_dat[[salinity]]
  warn_eos_range(t, s)

  density <- eos80_density(t, s)
  env_dat[[name %||% if (sigma) "sigma_theta" else "density"]] <-
    if (sigma) density - 1000 else density
  env_dat
}

#' UNESCO (1983) one-atmosphere equation of state
#'
#' Density of pure water from the Bigg (1967) polynomial, plus the salinity
#' terms. Written out in full rather than factored, so it can be checked line by
#' line against the published coefficients.
#'
#' @param t temperature in degrees C
#' @param s salinity in PSU
#' @return density in kg/m^3
#' @keywords internal
eos80_density <- function(t, s) {
  # Pure water.
  rho_w <- 999.842594 +
    6.793952e-2 * t -
    9.095290e-3 * t^2 +
    1.001685e-4 * t^3 -
    1.120083e-6 * t^4 +
    6.536332e-9 * t^5

  a <- 8.24493e-1 -
    4.0899e-3 * t +
    7.6438e-5 * t^2 -
    8.2467e-7 * t^3 +
    5.3875e-9 * t^4

  b <- -5.72466e-3 +
    1.0227e-4 * t -
    1.6546e-6 * t^2

  c <- 4.8314e-4

  rho_w + a * s + b * s^1.5 + c * s^2
}

#' Warn when temperature or salinity is outside the fitted range
#'
#' @param t temperature in degrees C
#' @param s salinity in PSU
#' @return invisible `NULL`
#' @keywords internal
warn_eos_range <- function(t, s) {
  # Salinity as a mass fraction rather than PSU sits comfortably inside the
  # fitted range, so the range check cannot see it. A whole column under 1 PSU
  # is not seawater, and the resulting density looks entirely plausible, which
  # is what makes it worth catching separately.
  usable <- s[!is.na(s)]
  if (length(usable) > 0 && max(usable) < 1) {
    warning(
      "Every salinity value is below 1 PSU, so the units are probably a mass ",
      "fraction rather than PSU.",
      "\n  Seawater near 0.035 would be taken as almost fresh, and the density ",
      "returned would be\n  about 1000 kg/m^3 rather than about 1027. Multiply ",
      "by 1000, or pass the PSU column.",
      call. = FALSE
    )
  }

  outside <- function(x, lo, hi) sum(!is.na(x) & (x < lo | x > hi))
  bad_t <- outside(t, -2, 40)
  bad_s <- outside(s, 0, 42)
  if (bad_t == 0 && bad_s == 0) return(invisible(NULL))

  parts <- c(
    if (bad_t > 0) paste0(bad_t, " temperature value(s) outside -2 to 40 C"),
    if (bad_s > 0) paste0(bad_s, " salinity value(s) outside 0 to 42 PSU")
  )
  warning(
    "The equation of state is being extrapolated: ",
    paste(parts, collapse = " and "), ".",
    "\n  Those rows get a number, but it is the value of a polynomial past the ",
    "data it was fitted to,\n  not a density. Check the units first: salinity ",
    "as a fraction rather than PSU is the\n  usual cause.",
    call. = FALSE
  )
  invisible(NULL)
}
