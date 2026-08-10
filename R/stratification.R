#' Buoyancy frequency between two depths
#'
#' How strongly the water column is stratified: the frequency at which a parcel
#' displaced vertically would oscillate, squared. Written \eqn{N^2}, and equal to
#' \eqn{-(g/\rho_0)\,\partial\rho/\partial z}.
#'
#' This is the quantity that decides how much work mixing has to do, and it
#' governs a great deal of what plankton do. A strongly stratified column keeps a
#' surface bloom in the light and cuts it off from nutrients below; a weakly
#' stratified one lets both mix. `vertical_gradient()` approximates the same idea
#' with a temperature difference, which is a reasonable proxy only where salinity
#' is uniform — and in the Gulf of Maine it is not, because Scotian Shelf inflow
#' is fresh enough to stratify water that is barely warmer at the surface.
#'
#' @section Getting the two levels:
#' `datamatch::accessEnvDat()` takes a `depth` argument, so a second call at a
#' different depth gives the deeper level, and the two are joined as columns on
#' the same points. Compute a density for each with [potential_density()], then
#' pass both here.
#'
#' ```r
#' surface <- datamatch::accessEnvDat(vars = c("SST", "SSS"), depth = c(0, 1), ...)
#' deep    <- datamatch::accessEnvDat(vars = c("SST", "SSS"), depth = c(90, 100), ...)
#' ```
#'
#' @section What a two-level estimate is and is not:
#' It is a **bulk** stratification over the layer between the two depths, not a
#' profile. Where a sharp pycnocline sits between the levels, the true peak
#' \eqn{N^2} at that interface is far larger than the layer average reported
#' here, and the depth of the peak is invisible.
#'
#' The choice of levels therefore does much of the work, and the trap is the
#' mixed layer. If the two straddle its base, \eqn{N^2} is dominated by how much
#' of the well-mixed layer was included rather than by the stratification
#' itself, and it will vary through the season for that reason alone. Compare
#' the levels against `MLD` before reading a seasonal cycle into the result.
#'
#' Negative values are returned rather than suppressed. They mean denser water
#' sits above lighter, which is genuine convective instability in winter and
#' otherwise usually a sign that the two levels are not what you think.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`
#' @param shallow name of the density column for the shallower level, in kg/m^3
#'   or as sigma-theta; both give the same answer since only the difference
#'   enters
#' @param deep name of the density column for the deeper level
#' @param depths the two depths in metres, shallower first, positive downwards
#' @param frequency return \eqn{N} in s^-1 instead of \eqn{N^2} in s^-2. `N` is
#'   `NA` wherever the column is unstable, since the root of a negative number
#'   is not a frequency
#' @param name name for the new column
#' @return `env_dat` with a stratification column
#' @examples
#' \dontrun{
#' env <- potential_density(env, "SST", "SSS", name = "rho_surface")
#' env <- potential_density(env, "SST_deep", "SSS_deep", name = "rho_deep")
#' env <- buoyancy_frequency(env, "rho_surface", "rho_deep",
#'                           depths = c(0.494, 92.326))
#' }
#' @seealso [potential_density()], [vertical_gradient()], [eady_growth_rate()]
#' @export
buoyancy_frequency <- function(env_dat, shallow, deep, depths,
                               frequency = FALSE, name = NULL) {
  resolve_vars(env_dat, c(shallow, deep))
  depths <- check_depths(depths)

  # Gravity, and a reference density for the Boussinesq approximation. Using a
  # constant rather than the local density changes the answer by well under a
  # percent and keeps the quantity comparable between places.
  gravity <- 9.81
  reference <- 1025

  # z is positive upwards and depth positive downwards, so the two sign flips
  # cancel: stable stratification, denser below, gives a positive N^2.
  spacing <- depths[2] - depths[1]
  n2 <- (gravity / reference) * (env_dat[[deep]] - env_dat[[shallow]]) / spacing

  warn_unstable(n2)
  env_dat[[name %||% if (frequency) "N" else "N2"]] <-
    if (frequency) sqrt(ifelse(n2 > 0, n2, NA_real_)) else n2
  env_dat
}

#' Eady growth rate of baroclinic instability
#'
#' How fast baroclinic instability grows: the rate at which the available
#' potential energy stored in tilted density surfaces is converted into eddies.
#' High where a sheared, weakly stratified flow can overturn, which on this
#' shelf means the shelf-break front and the edges of warm-core rings — the
#' places eddies are actually generated.
#'
#' Computed in the maximum-growth-rate form of Lindzen and Farrell (1980),
#' \eqn{\sigma = 0.31\,|f|\,|\partial U/\partial z| / N}, for the instability
#' described by Eady (1949).
#'
#' **Eady is a person, not a spelling of "eddy".** Eric Eady set out the model in
#' 1949. The collision is unlucky, because the rate named after him is precisely
#' a predictor of where eddies form, so the two words look interchangeable beside
#' each other and are not.
#'
#' @section What it needs:
#' Velocities at two depths and a stratification spanning the same layer, which
#' means two calls to `datamatch::accessEnvDat()` at different `depth` ranges
#' joined as columns. See [buoyancy_frequency()], which produces the
#' stratification and explains how to choose the levels.
#'
#' @section What it is and is not:
#' An index of **where and when** instability is favoured, not a prediction of
#' growth. The formula assumes quasi-geostrophic scaling and uniform
#' stratification through the layer, neither of which holds exactly on a shelf,
#' and the 0.31 is the maximum over wavenumber rather than the rate of any
#' particular disturbance. Read it as a comparative field.
#'
#' It is undefined where the column is not stably stratified, since \eqn{N} is
#' then not a frequency. Those cells come back `NA` rather than as a very large
#' rate, which is what dividing by a vanishing \eqn{N} would otherwise produce
#' and which would look like intense instability exactly where the assumption
#' has failed.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`
#' @param shallow length-2 character vector naming the eastward and northward
#'   velocity columns at the shallower level, in m/s
#' @param deep the same for the deeper level
#' @param depths the two depths in metres, shallower first
#' @param stratification name of an \eqn{N^2} column, in s^-2, as produced by
#'   [buoyancy_frequency()]
#' @param per `"day"` (the default) or `"second"`
#' @param name name for the new column
#' @return `env_dat` with a growth rate column, in day^-1 or s^-1
#' @references
#' Eady ET (1949). Long waves and cyclone waves. *Tellus* **1**(3), 33-52.
#' \doi{10.3402/tellusa.v1i3.8507}
#'
#' Lindzen RS, Farrell B (1980). A simple approximate result for the maximum
#' growth rate of baroclinic instabilities. *Journal of the Atmospheric
#' Sciences* **37**(7), 1648-1654.
#' \doi{10.1175/1520-0469(1980)037<1648:ASARFT>2.0.CO;2}
#' @examples
#' \dontrun{
#' env <- buoyancy_frequency(env, "rho_surface", "rho_deep",
#'                           depths = c(0.494, 92.326))
#' env <- eady_growth_rate(env,
#'                         shallow = c("UO", "VO"),
#'                         deep = c("UO_deep", "VO_deep"),
#'                         depths = c(0.494, 92.326))
#' }
#' @seealso [buoyancy_frequency()], [detect_eddies()], [flow_deformation()]
#' @export
eady_growth_rate <- function(env_dat, shallow = c("UO", "VO"), deep,
                             depths, stratification = "N2",
                             per = c("day", "second"), name = NULL) {
  per <- match.arg(per)
  if (length(shallow) != 2 || length(deep) != 2) {
    stop("`shallow` and `deep` must each name two columns, the eastward and ",
         "northward velocity at that level.", call. = FALSE)
  }
  resolve_vars(env_dat, c(shallow, deep, stratification))
  depths <- check_depths(depths)

  spacing <- depths[2] - depths[1]
  du <- (env_dat[[shallow[1]]] - env_dat[[deep[1]]]) / spacing
  dv <- (env_dat[[shallow[2]]] - env_dat[[deep[2]]]) / spacing
  shear <- sqrt(du^2 + dv^2)

  latitude <- sf::st_coordinates(env_dat)[, 2]
  f <- abs(coriolis_parameter(latitude))

  n2 <- env_dat[[stratification]]
  # N is only a frequency where the column is stably stratified.
  brunt <- sqrt(ifelse(n2 > 0, n2, NA_real_))

  rate <- 0.31 * f * shear / brunt
  if (per == "day") rate <- rate * 86400

  warn_eady_undefined(n2, f)
  env_dat[[name %||% "eady_growth"]] <- rate
  env_dat
}

#' Coriolis parameter from latitude
#'
#' @param latitude degrees north
#' @return `f` in s^-1, `NA` within two degrees of the equator where it vanishes
#'   and anything divided by it stops meaning anything
#' @keywords internal
coriolis_parameter <- function(latitude) {
  omega <- 7.2921e-5
  f <- 2 * omega * sin(latitude * pi / 180)
  ifelse(abs(latitude) < 2, NA_real_, f)
}

#' Validate a pair of depths
#'
#' @param depths two depths in metres
#' @return the depths, sorted shallow first
#' @keywords internal
check_depths <- function(depths) {
  if (!is.numeric(depths) || length(depths) != 2 || anyNA(depths)) {
    stop("`depths` must be two numbers, the depth of each level in metres.",
         call. = FALSE)
  }
  if (depths[1] == depths[2]) {
    stop("The two depths are the same, so there is no layer between them and ",
         "no vertical difference\n  to take.", call. = FALSE)
  }
  if (depths[1] > depths[2]) {
    stop("`depths` must be shallower first. Got ", depths[1], " then ",
         depths[2], ".\n  Swapping them silently would flip the sign of the ",
         "stratification.", call. = FALSE)
  }
  depths
}

#' Warn when the column is reported as unstable
#'
#' @param n2 buoyancy frequency squared
#' @return invisible `NULL`
#' @keywords internal
warn_unstable <- function(n2) {
  usable <- !is.na(n2)
  if (!any(usable)) return(invisible(NULL))
  share <- mean(n2[usable] < 0)
  if (share > 0.1) {
    warning(
      round(100 * share), "% of cells have denser water above lighter, so N^2 ",
      "is negative there.",
      "\n  That is real in winter convection, and otherwise usually means the ",
      "two levels are not\n  what they are taken to be -- check that the ",
      "shallower column really is the shallower\n  level, and that `depths` ",
      "matches the levels actually fetched.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Warn about cells where the growth rate cannot be defined
#'
#' @param n2 buoyancy frequency squared
#' @param f Coriolis parameter
#' @return invisible `NULL`
#' @keywords internal
warn_eady_undefined <- function(n2, f) {
  lost <- mean(is.na(n2) | n2 <= 0 | is.na(f))
  if (lost > 0.1) {
    warning(
      round(100 * lost), "% of cells have no growth rate, because the column ",
      "is not stably stratified\n  there or the latitude is too near the ",
      "equator for the Coriolis parameter to be usable.",
      "\n  Those are NA rather than a very large rate, which is what dividing ",
      "by a vanishing N would\n  give and which would look like intense ",
      "instability exactly where the assumption failed.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
