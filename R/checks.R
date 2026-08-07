#' Typical flow speed and domain size, for diagnosing empty Lagrangian output
#'
#' FTLE and FSLE both fail the same way: a domain that is small relative to how
#' far the flow carries a parcel leaves nothing to report. The numbers needed to
#' say so are the median speed and the extent, so they are gathered once here.
#'
#' @param env_dat an `sf` POINT object
#' @param u,v velocity column names, in m/s
#' @return list with `speed` (m/s), `width_km`, and `height_km`
#' @keywords internal
flow_scale <- function(env_dat, u, v) {
  speed <- sqrt(env_dat[[u]]^2 + env_dat[[v]]^2)
  xy <- sf::st_coordinates(env_dat)
  metres_per_degree <- 111320
  mean_lat <- mean(range(xy[, 2], na.rm = TRUE))

  list(
    speed = stats::median(speed, na.rm = TRUE),
    width_km = diff(range(xy[, 1], na.rm = TRUE)) * metres_per_degree *
      cos(mean_lat * pi / 180) / 1000,
    height_km = diff(range(xy[, 2], na.rm = TRUE)) * metres_per_degree / 1000
  )
}

#' Displacement a parcel undergoes over a window, in km
#'
#' @param speed speed in m/s
#' @param days window length in days
#' @return distance in km
#' @keywords internal
displacement_km <- function(speed, days) speed * 86400 * days / 1000

#' Warn when a Lagrangian diagnostic returns almost nothing
#'
#' A field that is nearly all `NA` is the expected result of asking for a longer
#' trajectory than the domain can hold, and every individual `NA` is correct.
#' What is not acceptable is returning that silently: the caller sees an empty
#' column and cannot tell an unsuitable domain from a broken function.
#'
#' The threshold is deliberately high. Losing a margin to the domain edge is
#' normal and does not need saying every time. Losing essentially everything
#' does.
#'
#' @param fn name of the calling function, for the message
#' @param na_fraction fraction of the output that is `NA`
#' @param detail one or more sentences diagnosing the cause
#' @param advice what the caller can change
#' @param threshold warn at or above this fraction
#' @return invisibly `NULL`
#' @keywords internal
warn_lagrangian <- function(fn, na_fraction, detail, advice, threshold = 0.9) {
  if (!is.finite(na_fraction) || na_fraction < threshold) return(invisible(NULL))

  headline <- if (na_fraction >= 1) {
    paste0(fn, "() returned no values at all: every point is NA.")
  } else {
    paste0(fn, "() returned almost nothing: ",
           format(round(100 * na_fraction, 1), nsmall = 1), "% of points are NA.")
  }

  warning(headline, "\n  ", detail, "\n  ", advice, call. = FALSE)
  invisible(NULL)
}

#' Is a column constant within every group?
#'
#' The shared test behind both degeneracy warnings: group by location and a
#' constant column is static through time, group by time step and it is uniform
#' across space.
#'
#' Implemented by counting distinct (group, value) pairs rather than by splitting
#' and looping, so it stays O(n) on the full point table. Comparison is exact:
#' two values differing in the last bit count as different, so this only fires on
#' columns that really are constant, never on a nearly-flat field.
#'
#' @param values a covariate column
#' @param group grouping vector of the same length
#' @return `TRUE` if every group holds a single distinct non-`NA` value
#' @keywords internal
constant_within <- function(values, group) {
  keep <- !is.na(values)
  if (!any(keep)) return(FALSE)
  pairs <- !duplicated(data.frame(group = group[keep], value = values[keep]))
  sum(pairs) == length(unique(group[keep]))
}

#' Warn about covariates a given operation cannot say anything about
#'
#' These are the cases that run to completion and return a column of perfectly
#' valid numbers carrying no information: the lag of a static covariate is the
#' covariate, the horizontal gradient of a basin-wide index is zero. Nothing
#' downstream can tell such a column from a real one, and a model handed it will
#' happily report it as uninformative rather than as a mistake — so the warning
#' is at the point where the intent is still visible.
#'
#' A warning rather than an error because the computation is well defined and the
#' caller may have meant it. `vars = NULL` expanding over an enriched object is
#' the case this is really aimed at.
#'
#' @param env_dat an `sf` POINT object
#' @param vars resolved covariate names
#' @param kind `"spatial"`, `"temporal"`, or `"any"` to skip the check
#' @return invisibly `NULL`; called for the warning
#' @keywords internal
warn_degenerate <- function(env_dat, vars, kind) {
  if (identical(kind, "any") || length(vars) == 0) return(invisible(NULL))

  if (identical(kind, "temporal")) {
    # With one time step every column is trivially constant through time, and a
    # temporal operation is already all-NA, so there is nothing to report.
    if (nrow(time_steps(env_dat)) < 2) return(invisible(NULL))

    location <- location_key(env_dat)
    flagged <- vars[vapply(vars, function(v) constant_within(env_dat[[v]], location),
                           logical(1))]
    if (length(flagged) > 0) {
      warning("Static covariate(s) in a temporal operation: ",
              paste(flagged, collapse = ", "), ".\n",
              "  These hold the same value at each location in every time step, ",
              "so a lag reproduces the column, a temporal gradient is zero, and ",
              "an integral is a running multiple of it.\n",
              "  Seafloor terrain from datamatch::attach_bathymetry() (DEPTH, ",
              "SLOPE, ASPECT, TPI) is static in this way. Name the variables you ",
              "meant, rather than letting `vars = NULL` take every column.",
              call. = FALSE)
    }
  }

  if (identical(kind, "spatial")) {
    step_index <- match_step_index(env_dat, time_steps(env_dat))
    flagged <- vars[vapply(vars, function(v) constant_within(env_dat[[v]], step_index),
                           logical(1))]
    if (length(flagged) > 0) {
      warning("Spatially uniform covariate(s) in a spatial operation: ",
              paste(flagged, collapse = ", "), ".\n",
              "  These take one value across the whole grid within each time ",
              "step, so a horizontal gradient is zero everywhere and a distance ",
              "to a contour of them is undefined.\n",
              "  Climate indices from datamatch::attach_climate_index() (NAO, ",
              "AO, AMO, PDO, LCR, AMOC) are uniform in this way: they carry ",
              "information about when, none about where.",
              call. = FALSE)
    }
  }

  invisible(NULL)
}
