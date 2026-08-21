#' Horizontal spatial gradient of a covariate
#'
#' Adds the magnitude of each covariate's horizontal gradient — how fast it
#' changes with distance — which is what picks out fronts, the convergence zones
#' where plankton aggregate.
#'
#' Computed by central differences on the covariate's own grid, in real distance
#' units rather than degrees. That distinction matters: a degree of longitude is
#' about 83 km at 42 degrees N but 111 km at the equator, so a gradient left in
#' per-degree units is stretched by latitude and not comparable across a study
#' area.
#'
#' This is a deliberate departure from `raster::terrain()`, which the original
#' pipeline of Ross et al. (2023) used for `sst_grad` and `uv_grad`.
#' `terrain()` treats its input as an
#' elevation in the same units as the coordinates and returns a slope angle, which
#' is dimensionally meaningless for a field measured in degrees Celsius. The
#' magnitude here is a real rate of change with a real unit.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param vars covariate columns to differentiate; `NULL` does all of them
#' @param per distance unit for the result: `"km"` (default) or `"m"`
#' @param components also add the eastward and northward components, as
#'   `<var>_grad_x` and `<var>_grad_y`
#' @param suffix suffix for the magnitude column
#' @return `env_dat` with a `<var>_grad` column per covariate, in covariate units
#'   per `per`
#' @references
#' Ross C, Runge J, Roberts J, Brady D, Tupper B, Record N (2023). Estimating
#' North Atlantic right whale prey based on Calanus finmarchicus thresholds.
#' *Marine Ecology Progress Series* **703**, 1-16. \doi{10.3354/meps14204}
#' @examples
#' \dontrun{
#' env <- datamatch::accessCopernicus(vars = c("SST", "SSS"), ...)
#' env <- horizontal_gradient(env, "SST")        # SST_grad, in degrees C per km
#' env <- horizontal_gradient(env, "SST", components = TRUE)
#'
#' # `vars = NULL` does every covariate column, which is right for a plain fetch
#' # and wrong once static or non-numeric columns have been attached.
#' env <- horizontal_gradient(env)
#' }
#' @export
horizontal_gradient <- function(env_dat, vars = NULL, per = c("km", "m"),
                                components = FALSE, suffix = "_grad") {
  per <- match.arg(per)
  vars <- resolve_vars(env_dat, vars, kind = "spatial")

  per_time_step(env_dat, vars, function(rast) {
    layers <- lapply(vars, function(v) {
      grad <- gradient_layers(rast[[v]], per = per)
      names(grad) <- paste0(v, c(suffix, paste0(suffix, "_x"), paste0(suffix, "_y")))
      if (components) grad else grad[[1]]
    })
    do.call(c, layers)
  })
}

#' Gradient magnitude and components of a single raster layer
#'
#' @param layer a one-layer `SpatRaster`
#' @param per distance unit, `"km"` or `"m"`
#' @return a three-layer `SpatRaster`: magnitude, eastward, northward
#' @keywords internal
gradient_layers <- function(layer, per = "km") {
  spacing <- cell_size(layer, per = per)

  # Central differences: (right - left) / 2dx and (up - down) / 2dy. The kernels
  # are written east-positive and north-positive, and terra's focal() runs rows
  # top-to-bottom (north first), so the y kernel is flipped relative to the x one.
  kernel_x <- matrix(c(0, 0, 0, -1, 0, 1, 0, 0, 0), nrow = 3, byrow = TRUE)
  kernel_y <- matrix(c(0, 1, 0, 0, 0, 0, 0, -1, 0), nrow = 3, byrow = TRUE)

  dz_dx <- terra::focal(layer, w = kernel_x, fun = sum, na.rm = FALSE) / (2 * spacing$x)
  dz_dy <- terra::focal(layer, w = kernel_y, fun = sum, na.rm = FALSE) / (2 * spacing$y)

  magnitude <- sqrt(dz_dx^2 + dz_dy^2)
  c(magnitude, dz_dx, dz_dy)
}

#' Cell size in real distance units
#'
#' Converts a lon/lat grid's degree spacing into distance. Longitude spacing
#' shrinks with the cosine of latitude, so it is computed per row rather than
#' once for the whole grid — at 45 degrees the error from ignoring that is about
#' 30 percent.
#'
#' A projected raster is already in linear units, so its resolution is used
#' directly.
#'
#' @param layer a `SpatRaster`
#' @param per `"km"` or `"m"`
#' @return `list(x =, y =)`, where `x` is a `SpatRaster` of per-row longitude
#'   spacing and `y` a scalar
#' @keywords internal
cell_size <- function(layer, per = "km") {
  scale <- if (per == "km") 1000 else 1
  resolution <- terra::res(layer)

  if (!terra::is.lonlat(layer)) {
    return(list(x = resolution[1] / scale, y = resolution[2] / scale))
  }

  # Mean Earth radius; good to a few tenths of a percent at any latitude.
  metres_per_degree <- 111320
  latitude <- terra::init(layer, "y")

  list(
    x = (resolution[1] * metres_per_degree * cos(latitude * pi / 180)) / scale,
    y = (resolution[2] * metres_per_degree) / scale
  )
}

#' Vertical temperature gradient
#'
#' The difference between surface and bottom temperature in each cell: a
#' stratification index, large where a warm surface layer sits over cold deep
#' water and near zero where the column is well mixed.
#'
#' Both defaults come from the same Copernicus physics dataset (`SST` and `BOTT`
#' in the `datamatch` catalog, `thetao` and `bottomT` as Copernicus codes), so
#' the usual case needs no extra fetch — it is a difference between two columns
#' already present.
#'
#' @section Other levels, and a better measure:
#' The two columns are arguments, not fixed, so this is not restricted to the
#' surface and the sea floor. `datamatch::accessCopernicus()` takes a `depth`
#' argument and returns one level per call, so a second fetch at a chosen depth
#' gives a temperature at that level, and the difference between any two levels
#' can be taken the same way.
#'
#' Where the two levels are known, [buoyancy_frequency()] is the better
#' quantity. This is a temperature difference, which stands in for
#' stratification only where salinity is uniform — and in the Gulf of Maine it
#' is not, since Scotian Shelf inflow is fresh enough to stratify water barely
#' warmer at the surface. \eqn{N^2} counts both contributions with the weights
#' the equation of state gives them.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param surface name of the shallower temperature column
#' @param bottom name of the deeper temperature column
#' @param depth optional depth column name; when given, the difference is divided
#'   by depth to give a per-metre rate rather than a total difference.
#'   `datamatch::attach_bathymetry()` supplies a `DEPTH` column.
#' @param name name for the new column
#' @return `env_dat` with the vertical gradient column added
#' @examples
#' \dontrun{
#' env <- datamatch::accessCopernicus(vars = c("SST", "BOTT"), ...)
#' env <- vertical_gradient(env)                 # SST and BOTT are the defaults
#'
#' # As a per-metre rate, with depth from datamatch's bathymetry
#' bathy <- datamatch::fetch_bathymetry(bounding_box = bb)
#' env <- datamatch::attach_bathymetry(env, bathy, "DEPTH")
#' env <- vertical_gradient(env, depth = "DEPTH")
#'
#' # Between two chosen levels rather than surface and sea floor
#' deep <- datamatch::accessCopernicus(vars = "SST", depth = c(90, 100), ...)
#' env$SST_90m <- deep$SST
#' env <- vertical_gradient(env, surface = "SST", bottom = "SST_90m")
#' }
#' @seealso [buoyancy_frequency()], [potential_density()]
#' @export
vertical_gradient <- function(env_dat, surface = "SST", bottom = "BOTT",
                              depth = NULL, name = NULL) {
  resolve_vars(env_dat, c(surface, bottom, depth))

  difference <- env_dat[[surface]] - env_dat[[bottom]]

  if (!is.null(depth)) {
    # Guard against dividing by a zero or negative depth, which would produce
    # Inf or a sign flip rather than a missing value.
    divisor <- env_dat[[depth]]
    divisor[!is.na(divisor) & divisor <= 0] <- NA_real_
    difference <- difference / divisor
  }

  env_dat[[name %||% paste0(surface, "_", bottom, "_vgrad")]] <- difference
  env_dat
}

#' Temporal gradient of a covariate
#'
#' Rate of change between consecutive time steps at each location — how fast
#' conditions are shifting, as distinct from what they currently are. A water mass
#' warming rapidly is a different habitat from one sitting at the same temperature.
#'
#' The first time step has no predecessor and is `NA`.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param vars covariate columns to differentiate; `NULL` does all of them
#' @param per time unit for the rate: `"step"` (default, change per time step),
#'   `"day"`, or `"month"`
#' @param suffix suffix for the new columns
#' @return `env_dat` with a `<var>_tgrad` column per covariate
#' @examples
#' \dontrun{
#' env <- temporal_gradient(env, "SST")                 # change per time step
#' env <- temporal_gradient(env, "SST", per = "day")    # degrees C per day
#' }
#' @export
temporal_gradient <- function(env_dat, vars = NULL,
                              per = c("step", "day", "month"), suffix = "_tgrad") {
  per <- match.arg(per)
  # No `kind` here on purpose: the lag_covariate() call below resolves the same
  # vars as "temporal" and issues the one warning. Setting it here too would
  # warn twice for a single user-facing call.
  vars <- resolve_vars(env_dat, vars)

  lagged <- lag_covariate(env_dat, vars, n = 1, suffix = "__prev")
  steps <- time_steps(env_dat)
  elapsed <- step_spacing(steps, per = per)

  for (v in vars) {
    change <- lagged[[v]] - lagged[[paste0(v, "__prev")]]
    divisor <- elapsed[match_step_index(env_dat, steps)]
    env_dat[[paste0(v, suffix)]] <- change / divisor
  }
  env_dat
}

#' Elapsed time between consecutive steps
#'
#' @param steps a time-step table from `time_steps()`
#' @param per `"step"`, `"day"`, or `"month"`
#' @return numeric vector, one per step; the first is `NA`
#' @keywords internal
step_spacing <- function(steps, per = "step") {
  if (per == "step") return(c(NA_real_, rep(1, nrow(steps) - 1)))

  dates <- as.Date(paste(steps$YEAR, steps$MONTH, steps$DAY, sep = "-"))
  days <- c(NA_real_, as.numeric(diff(dates)))
  if (per == "day") days else days / 30.4375  # mean month length
}

#' Index of each row's time step within the ordered step table
#'
#' @param env_dat an `sf` POINT object
#' @param steps a time-step table from `time_steps()`
#' @return integer vector, one per row of `env_dat`
#' @keywords internal
match_step_index <- function(env_dat, steps) {
  key <- function(y, m, d) paste(y, m, d, sep = "-")
  match(key(env_dat$YEAR, env_dat$MONTH, env_dat$DAY),
        key(steps$YEAR, steps$MONTH, steps$DAY))
}
