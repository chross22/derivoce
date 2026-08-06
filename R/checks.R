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
              "AO, AMO, PDO) are uniform in this way: they carry information ",
              "about when, none about where.",
              call. = FALSE)
    }
  }

  invisible(NULL)
}
