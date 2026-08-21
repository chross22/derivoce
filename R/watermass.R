#' Water-mass fraction from two endmembers
#'
#' Estimates what fraction of each cell is the first of two named water masses,
#' by projecting its temperature and salinity onto the mixing line between the
#' two endmembers in T-S space.
#'
#' The Gulf of Maine is the standard setting for this. Its deep water is a
#' mixture of Warm Slope Water and Labrador Slope Water in proportions that vary
#' episodically, and the proportion sets the nutrient supply: Labrador Slope
#' Water is colder, fresher, and nutrient-poor, Warm Slope Water warmer, saltier,
#' and nutrient-rich (Townsend et al. 2015). Surface water is instead a mixture
#' of Scotian Shelf Water with slope water, which is the same calculation with
#' different endmembers.
#'
#' @section How the fraction is computed:
#' A cell's `(T, S)` is projected onto the line joining the two endmembers, and
#' the fraction is where it lands: 1 at the first endmember, 0 at the second.
#'
#' Temperature and salinity are on different scales, so each axis is normalised
#' by the endmember separation in that variable before projecting. Without that
#' the axis with the larger numerical spread would dominate the answer for
#' reasons of units rather than of oceanography.
#'
#' Values are clamped to `[0, 1]`. A cell beyond an endmember is either a third
#' water mass or an endmember choice that does not bracket the data, and reporting
#' a fraction of 1.4 would dress up that problem as a measurement.
#'
#' @section Check the residual:
#' Projection always returns a fraction, even for water that is not a mixture of
#' these two masses at all. `residual = TRUE` adds the distance from the mixing
#' line in the same normalised units, which is the number that says whether the
#' fraction means anything: near zero, the cell sits on the line and the mixture
#' is a good description; large, and it does not.
#'
#' This is why the endmembers are an argument with no default. They vary by
#' region, season, and year, and a wrong pair produces confident nonsense.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param endmembers a list of exactly two named `c(temperature, salinity)`
#'   vectors. The fraction returned is of the **first**.
#' @param temperature name of the temperature column, in degrees C
#' @param salinity name of the salinity column, in PSU
#' @param residual also add the distance from the mixing line, as
#'   `<name>_residual`
#' @param name name for the fraction column. `NULL` uses the first endmember's
#'   name, lowercased, with `_frac`.
#' @return `env_dat` with a fraction column in `[0, 1]`, and optionally a
#'   residual column
#' @references
#' Townsend DW, Pettigrew NR, Thomas MA, Neary MG, McGillicuddy DJ, O'Donnell J
#' (2015). Water masses and nutrient sources to the Gulf of Maine. *Journal of
#' Marine Research* **73**, 93-122.
#' @examples
#' \dontrun{
#' # Deep water: which slope water is filling the basins?
#' env <- water_mass_fraction(env, endmembers = list(
#'   LSW = c(temperature = 6,  salinity = 34.4),
#'   WSW = c(temperature = 12, salinity = 35.4)
#' ), residual = TRUE)
#'
#' # Surface water: how much of this is Scotian Shelf Water?
#' env <- water_mass_fraction(env, endmembers = list(
#'   SSW   = c(temperature = 2,  salinity = 32.0),
#'   slope = c(temperature = 10, salinity = 35.0)
#' ))
#' }
#' @seealso [derived_indices()]
#' @export
water_mass_fraction <- function(env_dat, endmembers, temperature = "SST",
                                salinity = "SSS", residual = FALSE,
                                name = NULL) {
  resolve_vars(env_dat, c(temperature, salinity))
  endmembers <- check_endmembers(endmembers)

  first <- endmembers[[1]]
  second <- endmembers[[2]]

  # Normalise each axis by the endmember separation, so temperature and salinity
  # contribute on comparable footing rather than in proportion to their units.
  scale_t <- abs(first[1] - second[1])
  scale_s <- abs(first[2] - second[2])

  observed_t <- (env_dat[[temperature]] - second[1]) / scale_t
  observed_s <- (env_dat[[salinity]] - second[2]) / scale_s
  axis_t <- (first[1] - second[1]) / scale_t
  axis_s <- (first[2] - second[2]) / scale_s

  projection <- (observed_t * axis_t + observed_s * axis_s) /
    (axis_t^2 + axis_s^2)

  column <- name %||% paste0(tolower(names(endmembers)[1]), "_frac")
  env_dat[[column]] <- pmin(pmax(projection, 0), 1)

  if (residual) {
    # Perpendicular distance from the mixing line, before clamping. Water that
    # is not a mixture of these two masses still gets a fraction; this is what
    # reveals it.
    env_dat[[paste0(column, "_residual")]] <- sqrt(
      (observed_t - projection * axis_t)^2 +
        (observed_s - projection * axis_s)^2
    )
  }

  env_dat
}

#' Validate a pair of T-S endmembers
#'
#' @param endmembers a candidate list of two named numeric pairs
#' @return the endmembers, with names filled in if absent
#' @keywords internal
check_endmembers <- function(endmembers) {
  if (!is.list(endmembers) || length(endmembers) != 2) {
    stop("`endmembers` must be a list of exactly two c(temperature, salinity) ",
         "vectors. Mixing more than two water masses needs a different method ",
         "than projection onto a line.", call. = FALSE)
  }

  ok <- vapply(endmembers, function(e) {
    is.numeric(e) && length(e) == 2 && all(is.finite(e))
  }, logical(1))
  if (!all(ok)) {
    stop("Each endmember must be a numeric c(temperature, salinity).",
         call. = FALSE)
  }

  if (is.null(names(endmembers)) || any(!nzchar(names(endmembers)))) {
    names(endmembers) <- c("first", "second")
  }

  first <- endmembers[[1]]
  second <- endmembers[[2]]
  if (isTRUE(all.equal(first[1], second[1])) &&
      isTRUE(all.equal(first[2], second[2]))) {
    stop("The two endmembers are identical, so there is no mixing line.",
         call. = FALSE)
  }
  if (isTRUE(all.equal(first[1], second[1]))) {
    stop("The endmembers have the same temperature, so temperature carries no ",
         "information here. Use a salinity threshold instead, or pick ",
         "endmembers that differ in both.", call. = FALSE)
  }
  if (isTRUE(all.equal(first[2], second[2]))) {
    stop("The endmembers have the same salinity, so salinity carries no ",
         "information here. Use a temperature threshold instead, or pick ",
         "endmembers that differ in both.", call. = FALSE)
  }

  endmembers
}
