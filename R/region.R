#' Anomaly of a covariate averaged over a box
#'
#' Averages a covariate over a lon/lat box in each time step, then subtracts a
#' reference to give an anomaly. The result is **one number per time step**,
#' broadcast to every row, so it behaves like a climate index rather than a map.
#'
#' This is the simplest of the region-scale indices and often the most robust.
#' It asks whether conditions in a place were unusual, without needing velocities
#' or endmembers, so it survives on products where the other methods cannot run.
#'
#' @section What it cannot tell you:
#' A box mean says conditions changed, not that water moved. A fresh anomaly in
#' the eastern Gulf of Maine is consistent with more Scotian Shelf inflow, and
#' also with local runoff, rainfall, or ice melt. Where a transport across a
#' section measures the crossing directly, this measures its most visible
#' consequence and asks you to supply the interpretation.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`
#' @param var covariate column to average
#' @param box named list with `xmin`, `xmax`, `ymin`, `ymax`, in degrees
#' @param reference `"climatology"` (the default) removes a separate mean per
#'   calendar month, so only departures from the usual conditions for that month
#'   survive. `"record"` removes one mean over the whole series, leaving the
#'   seasonal cycle in. `"none"` returns the box mean itself.
#' @param name name for the new column
#' @return `env_dat` with an anomaly column added, in the units of `var`
#' @examples
#' \dontrun{
#' env <- box_anomaly(env, "SSS", box = list(
#'   xmin = -68.5, xmax = -66.5, ymin = 43.0, ymax = 44.5
#' ))
#' }
#' @seealso [eastern_gom_salinity()], [derived_indices()]
#' @export
box_anomaly <- function(env_dat, var, box,
                        reference = c("climatology", "record", "none"),
                        name = NULL) {
  reference <- match.arg(reference)
  resolve_vars(env_dat, var)
  box <- check_box(box)

  xy <- sf::st_coordinates(env_dat)
  inside <- xy[, 1] >= box$xmin & xy[, 1] <= box$xmax &
    xy[, 2] >= box$ymin & xy[, 2] <= box$ymax

  if (!any(inside)) {
    stop("The box contains no grid points.",
         "\n  box:  x [", box$xmin, ", ", box$xmax, "]  y [",
         box$ymin, ", ", box$ymax, "]",
         "\n  data: x [", signif(min(xy[, 1]), 6), ", ", signif(max(xy[, 1]), 6),
         "]  y [", signif(min(xy[, 2]), 6), ", ", signif(max(xy[, 2]), 6), "]",
         call. = FALSE)
  }

  steps <- time_steps(env_dat)
  step_index <- match_step_index(env_dat, steps)
  values <- env_dat[[var]]

  # One mean per time step, over the cells inside the box.
  box_mean <- vapply(seq_len(nrow(steps)), function(i) {
    mean(values[step_index == i & inside], na.rm = TRUE)
  }, numeric(1))

  anomaly <- switch(
    reference,
    none = box_mean,
    record = box_mean - mean(box_mean, na.rm = TRUE),
    climatology = box_mean - stats::ave(box_mean, steps$MONTH,
                                        FUN = function(z) mean(z, na.rm = TRUE))
  )

  suffix <- switch(reference, none = "_box", record = "_box_anom",
                   climatology = "_box_anom")
  env_dat[[name %||% paste0(var, suffix)]] <- anomaly[step_index]
  env_dat
}

#' Eastern Gulf of Maine salinity index
#'
#' Surface salinity anomaly over the eastern Gulf of Maine, the region where
#' Scotian Shelf inflow first appears after rounding Cape Sable. Negative values
#' are fresher than usual, which indicates a stronger inflow of cold fresh
#' Scotian Shelf Water.
#'
#' Follows the approach of Grodsky et al. (2025), who showed that satellite
#' surface salinity in the eastern Gulf tracks winter Scotian Shelf inflow and
#' relates it to the coastal and interior pathways the water then takes. Their
#' index is built from SMAP satellite salinity; this computes the same kind of
#' quantity from whatever gridded salinity field you supply, so it is **not a
#' reproduction of their published series** and should not be compared to it
#' value for value.
#'
#' @section Why this is the seasonal one:
#' The signal is a winter one. Scotian Shelf Water arrives in the upper layers of
#' the Gulf in winter, roughly six months after the late-summer deep influx of
#' warm salty water. A whole-year mean mixes the two regimes together and dilutes
#' the thing being measured, so restrict to winter months at the fetch, or subset
#' afterwards, rather than reading an annual value as an inflow index.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`
#' @param salinity name of the surface salinity column, in PSU
#' @param reference passed to [box_anomaly()]
#' @param name name for the new column
#' @return `env_dat` with an `egom_salinity` column, in PSU
#' @references
#' Grodsky SA, Vandemark D, Levin J (2025). An eastern Gulf of Maine salinity
#' index for monitoring winter Scotian Shelf inflow and its relation to coastal
#' and interior pathways. *Journal of Geophysical Research: Oceans* **130**(5).
#' \doi{10.1029/2024JC021891}
#' @examples
#' \dontrun{
#' env <- eastern_gom_salinity(env)
#' }
#' @seealso [box_anomaly()], [eastern_gom_box()], [scotian_shelf_inflow()]
#' @export
eastern_gom_salinity <- function(env_dat, salinity = "SSS",
                                 reference = c("climatology", "record", "none"),
                                 name = "egom_salinity") {
  box_anomaly(env_dat, var = salinity, box = eastern_gom_box(),
              reference = match.arg(reference), name = name)
}

#' The eastern Gulf of Maine box
#'
#' Exposed so the region behind [eastern_gom_salinity()] can be inspected or
#' plotted. Approximate, and not the exact footprint of any published index.
#'
#' @return a named list with `xmin`, `xmax`, `ymin`, `ymax`
#' @examples
#' eastern_gom_box()
#' @export
eastern_gom_box <- function() {
  list(xmin = -68.0, xmax = -66.0, ymin = 43.0, ymax = 44.5)
}

#' Validate a bounding box
#'
#' @param box a candidate named list
#' @return the box, as a plain list
#' @keywords internal
check_box <- function(box) {
  needed <- c("xmin", "xmax", "ymin", "ymax")
  if (!is.list(box) || !all(needed %in% names(box))) {
    stop("`box` must be a list with xmin, xmax, ymin, and ymax.", call. = FALSE)
  }
  box <- lapply(box[needed], as.numeric)
  if (any(!is.finite(unlist(box)))) {
    stop("`box` values must all be finite numbers.", call. = FALSE)
  }
  if (box$xmin >= box$xmax || box$ymin >= box$ymax) {
    stop("`box` is empty: xmin must be below xmax and ymin below ymax.",
         call. = FALSE)
  }
  box
}
