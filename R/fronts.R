#' Distance to the nearest front
#'
#' Fronts are where a covariate changes sharply over a short distance, and they
#' concentrate plankton. But the local gradient alone is a blunt predictor of
#' that: a station sitting in the middle of a smooth patch has a gradient of zero
#' whether the nearest front is 2 km away or 200 km away, and those are very
#' different places to be. Distance to the nearest front separates them.
#'
#' Fronts are identified by thresholding the covariate's horizontal gradient, then
#' the distance from every cell to the nearest front cell is computed with a
#' distance transform.
#'
#' @section Choosing the threshold:
#' `quantile` (the default route) calls a cell frontal if its gradient is in the
#' top fraction of the distribution — `0.9` means the strongest 10% of gradients.
#'
#' `scope` decides what that fraction is taken over, and the two answers mean
#' different things:
#'
#' \itemize{
#'   \item `"record"` (the default) takes one threshold across the whole time
#'     series, so "front" means the same physical sharpness in every month. Winter
#'     months with weak stratification may then contain no fronts at all — which is
#'     a real result, not a failure.
#'   \item `"step"` takes a separate threshold per time step, so every month has
#'     fronts by construction. Use this when the question is where the sharpest
#'     features are *this month*, and not whether this month has sharp features.
#' }
#'
#' Pass `threshold` instead to set an absolute gradient value, which is the right
#' choice when a physically meaningful cutoff is known.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param var covariate whose gradient defines the fronts
#' @param threshold absolute gradient value above which a cell is frontal; takes
#'   precedence over `quantile` when given
#' @param quantile fraction of the gradient distribution treated as frontal
#' @param scope `"record"` or `"step"`; see above
#' @param per distance unit for the result: `"km"` (default) or `"m"`
#' @param name name for the new column
#' @return `env_dat` with a distance-to-front column. Cells on a front are `0`. A
#'   time step containing no front is `NA` throughout, since "distance to nothing"
#'   has no value.
#' @references
#' Identifying fronts by thresholding a gradient field follows the standard
#' approach, of which Belkin and O'Reilly's is the reference implementation for
#' satellite SST and chlorophyll. Theirs adds a contextual median filter that
#' preserves front shape before the gradient is taken, which matters on noisy
#' satellite retrievals and less on a model field that is already smooth. What is
#' done here is the plain gradient threshold, with `scope` deciding whether the
#' threshold is one value for the record or one per time step.
#'
#' Belkin IM, O'Reilly JE (2009). An algorithm for oceanic front detection in
#' chlorophyll and SST satellite imagery. *Journal of Marine Systems* **78**(3),
#' 319-326. \doi{10.1016/j.jmarsys.2008.11.018}
#' @examples
#' \dontrun{
#' # Distance to the sharpest 10% of thermal gradients
#' env <- distance_to_front(env, "SST")
#' # A physically defined front: 0.05 degrees C per km
#' env <- distance_to_front(env, "SST", threshold = 0.05)
#' }
#' @export
distance_to_front <- function(env_dat, var, threshold = NULL, quantile = 0.9,
                              scope = c("record", "step"), per = c("km", "m"),
                              name = NULL) {
  scope <- match.arg(scope)
  per <- match.arg(per)
  # No `kind` here on purpose: the horizontal_gradient() call below resolves the
  # same var as "spatial" and issues the one warning.
  resolve_vars(env_dat, var)

  if (!is.null(threshold) && threshold <= 0) {
    stop("threshold must be positive; a gradient magnitude is never negative.",
         call. = FALSE)
  }
  if (is.null(threshold) && (quantile <= 0 || quantile >= 1)) {
    stop("quantile must be strictly between 0 and 1.", call. = FALSE)
  }

  gradient_column <- paste0(var, "__front_grad")
  with_gradient <- horizontal_gradient(env_dat, var, per = per,
                                       suffix = "__front_grad")
  gradients <- with_gradient[[gradient_column]]

  record_threshold <- threshold %||% stats::quantile(gradients, quantile,
                                                      na.rm = TRUE, names = FALSE)

  steps <- time_steps(env_dat)
  result <- rep(NA_real_, nrow(env_dat))
  scale <- if (per == "km") 1000 else 1

  for (i in seq_len(nrow(steps))) {
    rows <- step_rows(env_dat, steps[i, ])
    step_gradients <- gradients[rows]

    cutoff <- if (!is.null(threshold) || scope == "record") {
      record_threshold
    } else {
      stats::quantile(step_gradients, quantile, na.rm = TRUE, names = FALSE)
    }

    result[rows] <- step_front_distance(with_gradient[rows, ], gradient_column,
                                        cutoff, scale)
  }

  env_dat[[name %||% paste0(var, "_front_dist")]] <- result
  env_dat
}

#' Distance to a contour of a covariate
#'
#' Distance from every point to the nearest place where a covariate crosses a
#' given value. Position relative to a contour is often what matters rather than
#' the value itself — plankton distributions track the shelf break, and "20 km
#' inshore of the 100 m isobath" locates that far better than "depth = 85 m" does.
#'
#' A cell is treated as on the contour if `level` falls between its value and
#' that of any neighbour, so the contour is found even though no cell sits
#' exactly on it.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param var covariate whose contours are wanted
#' @param levels one or more values to contour
#' @param per distance unit for the result: `"km"` (default) or `"m"`
#' @param prefix stem for the new column names; defaults to `var`
#' @return `env_dat` with one distance column per level. A level lying entirely
#'   outside the covariate's range in a time step gives `NA` for that step.
#' @examples
#' \dontrun{
#' # Distance to the shelf break and two other isobaths
#' env <- distance_to_contour(env, "DEPTH", levels = c(50, 100, 200))
#' # -> DEPTH_dist_50, DEPTH_dist_100, DEPTH_dist_200
#'
#' # Works on any field: distance to the 10 degree isotherm
#' env <- distance_to_contour(env, "SST", levels = 10)
#' }
#' @export
distance_to_contour <- function(env_dat, var, levels, per = c("km", "m"),
                                prefix = NULL) {
  per <- match.arg(per)
  resolve_vars(env_dat, var, kind = "spatial")

  if (length(levels) == 0 || any(!is.finite(levels))) {
    stop("levels must be one or more finite values.", call. = FALSE)
  }

  steps <- time_steps(env_dat)
  scale <- if (per == "km") 1000 else 1
  prefix <- prefix %||% var

  for (level in levels) {
    result <- rep(NA_real_, nrow(env_dat))

    for (i in seq_len(nrow(steps))) {
      rows <- step_rows(env_dat, steps[i, ])
      points <- env_dat[rows, ]
      values <- points[[var]]

      # The contour has to be bracketed by the data; a level outside the range
      # entirely has no location in this time step.
      if (all(is.na(values)) || level < min(values, na.rm = TRUE) ||
          level > max(values, na.rm = TRUE)) {
        next
      }

      grid <- rasterize_step(points, var)
      on_contour <- contour_cells(grid, level)
      if (!any(terra::values(on_contour) == 1, na.rm = TRUE)) next

      distances <- terra::distance(on_contour)
      result[rows] <- terra::extract(distances, sf::st_coordinates(points))[, 1] / scale
    }

    env_dat[[paste0(prefix, "_dist_", level)]] <- result
  }
  env_dat
}

#' Distance to isobaths
#'
#' Convenience wrapper on [distance_to_contour()] for depth contours, the most
#' common use.
#'
#' @param env_dat an `sf` POINT object carrying a depth column.
#'   `datamatch::attach_bathymetry()` adds one, named `DEPTH`.
#' @param depth name of the depth column, as a positive magnitude
#' @param levels isobaths to measure to, in the units of `depth`
#' @param per distance unit for the result
#' @return `env_dat` with one distance column per isobath
#' @examples
#' \dontrun{
#' bathy <- datamatch::fetch_bathymetry(bounding_box = bb)
#' env <- datamatch::attach_bathymetry(env, bathy, "DEPTH")
#' env <- distance_to_isobath(env, levels = c(50, 100, 200))
#' }
#' @export
distance_to_isobath <- function(env_dat, depth = "DEPTH",
                                levels = c(50, 100, 200), per = c("km", "m")) {
  distance_to_contour(env_dat, depth, levels, per = match.arg(per),
                      prefix = "isobath")
}

#' Cells a contour passes through
#'
#' A contour almost never falls exactly on a cell centre, so a cell is marked
#' when the level lies between its own value and any of its neighbours' — that
#' is, when the field crosses the level somewhere in the cell's neighbourhood.
#'
#' @param grid a one-layer `SpatRaster`
#' @param level the value being contoured
#' @return a `SpatRaster` of `1` on the contour and `NA` elsewhere, as
#'   `terra::distance()` expects
#' @keywords internal
contour_cells <- function(grid, level) {
  above <- grid >= level
  # A 3x3 range of the above/below indicator is 1 exactly where the
  # neighbourhood contains both sides of the level, i.e. where it is crossed.
  lowest <- terra::focal(above, w = 3, fun = min, na.rm = TRUE)
  highest <- terra::focal(above, w = 3, fun = max, na.rm = TRUE)

  terra::ifel(highest > lowest, 1, NA)
}

#' Distance to the nearest frontal cell within one time step
#'
#' @param points an `sf` POINT object for a single time step, carrying a
#'   gradient column
#' @param gradient_column name of the gradient column
#' @param cutoff gradient value at or above which a cell is frontal
#' @param scale divisor converting metres to the requested unit
#' @return numeric vector of distances, one per point
#' @keywords internal
step_front_distance <- function(points, gradient_column, cutoff, scale) {
  gradients <- points[[gradient_column]]
  if (all(is.na(gradients)) || !any(gradients >= cutoff, na.rm = TRUE)) {
    return(rep(NA_real_, nrow(points)))
  }

  grid <- rasterize_step(points, gradient_column)
  # terra::distance() measures from NA cells to the nearest non-NA cell, so the
  # fronts are what must be left non-NA.
  fronts <- terra::ifel(grid >= cutoff, 1, NA)

  distances <- terra::distance(fronts)
  terra::extract(distances, sf::st_coordinates(points))[, 1] / scale
}
