#' Finite-Time Lyapunov Exponent
#'
#' Measures how fast two initially neighbouring water parcels separate over an
#' integration window. Ridges of high FTLE mark Lagrangian coherent structures —
#' the transport barriers and convergence lines that organise surface flow.
#'
#' The direction matters, and is the reason this takes an argument rather than
#' picking one:
#'
#' \itemize{
#'   \item **Backward** (the default) integrates into the past, so ridges are
#'     *attracting* structures: places where water arriving from different origins
#'     converges. These are where buoyant particles and weak swimmers pile up, so
#'     this is usually the direction of interest for plankton aggregation.
#'   \item **Forward** integrates into the future, so ridges are *repelling*
#'     structures: barriers that neighbouring parcels straddle and are pulled
#'     apart by. These separate water masses rather than concentrating them.
#' }
#'
#' Particles are seeded on the covariate grid at each time step, advected through
#' the time-varying velocity field with a fourth-order Runge-Kutta scheme, and the
#' resulting flow map is differentiated to give the Cauchy-Green strain tensor.
#' The FTLE is then
#'
#' \deqn{\sigma = \frac{1}{|T|} \ln \sqrt{\lambda_{max}(C)}}
#'
#' where \eqn{\lambda_{max}(C)} is the largest eigenvalue of the strain tensor and
#' \eqn{T} the integration time. Units are inverse days.
#'
#' @section Choosing an integration time:
#' `integration_days` sets the scale of structure resolved: too short and no
#' separation has accumulated, too long and the ridges blur as particles wander
#' far from where they started. For mesoscale ocean features, days to a few weeks
#' is the usual range.
#'
#' Note that monthly-mean velocity fields, which is what
#' `datamatch::accessEnvDat()` returns for a `P1M` product, have already averaged
#' away the eddies that generate the sharpest structures. FTLE computed from
#' monthly means is a real diagnostic of the mean circulation, but it is not the
#' same quantity as FTLE from daily or hourly fields, and its ridges are
#' correspondingly smoother. Prefer a daily (`P1D`) product where the choice
#' exists.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`, on a
#'   regular lon/lat grid, containing eastward and northward velocity columns
#' @param u name of the eastward velocity column, in m/s
#' @param v name of the northward velocity column, in m/s
#' @param integration_days length of the integration window, in days
#' @param direction `"backward"` (attracting structures, the default) or
#'   `"forward"` (repelling structures)
#' @param step_hours Runge-Kutta step size, in hours. Smaller is more accurate and
#'   slower; the default resolves a particle crossing roughly one grid cell per
#'   step at typical shelf speeds.
#' @param name name for the new column
#' @return `env_dat` with an FTLE column added, in units of 1/day. Grid-boundary
#'   cells and particles that leave the domain or run aground are `NA`.
#' @examples
#' \dontrun{
#' env <- datamatch::accessEnvDat(
#'   dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
#'   vars = c("UO", "VO"), ...
#' )
#' # Where water converges - candidate aggregation sites
#' env <- ftle(env, integration_days = 14)       # UO and VO are the defaults
#' # Where water masses are pulled apart - transport barriers
#' env <- ftle(env, direction = "forward")
#'
#' # At depth: one model level per fetch
#' deep <- datamatch::accessEnvDat(vars = c("UO", "VO"), depth = c(100, 100), ...)
#' deep <- ftle(deep, integration_days = 14)
#' }
#' @references
#' Haller, G. (2015). Lagrangian coherent structures. *Annual Review of Fluid
#' Mechanics*, 47, 137-162.
#' @export
ftle <- function(env_dat, u = "UO", v = "VO", integration_days = 14,
                 direction = c("backward", "forward"), step_hours = 6,
                 name = NULL) {
  direction <- match.arg(direction)
  resolve_vars(env_dat, c(u, v))

  if (integration_days <= 0) {
    stop("integration_days must be positive; use `direction` to integrate ",
         "into the past.", call. = FALSE)
  }
  if (step_hours <= 0 || step_hours > integration_days * 24) {
    stop("step_hours must be positive and no longer than the integration window.",
         call. = FALSE)
  }

  steps <- time_steps(env_dat)
  step_times <- as.numeric(as.Date(paste(steps$YEAR, steps$MONTH, steps$DAY,
                                          sep = "-")))
  velocity <- lapply(seq_len(nrow(steps)), function(i) {
    rasterize_step(env_dat[step_rows(env_dat, steps[i, ]), ], c(u, v))
  })

  sign <- if (direction == "backward") -1 else 1
  result <- rep(NA_real_, nrow(env_dat))
  left_domain <- 0
  particles <- 0

  for (i in seq_len(nrow(steps))) {
    rows <- step_rows(env_dat, steps[i, ])
    seeds <- sf::st_coordinates(env_dat[rows, ])

    final <- advect(seeds, start_time = step_times[i],
                    duration = sign * integration_days,
                    step_days = step_hours / 24,
                    velocity = velocity, times = step_times)

    # A particle whose trajectory ran outside the velocity field comes back NA,
    # and poisons the central differences of its neighbours too. Counting them
    # is what lets an empty result be explained rather than just returned.
    left_domain <- left_domain + sum(!stats::complete.cases(final))
    particles <- particles + nrow(final)

    result[rows] <- ftle_from_flow_map(seeds, final, integration_days)
  }

  scale <- flow_scale(env_dat, u, v)
  reach <- displacement_km(scale$speed, integration_days)
  edge <- if (direction == "backward") "upstream" else "downstream"

  warn_lagrangian(
    "ftle", mean(is.na(result)),
    detail = paste0(
      "A ", integration_days, "-day integration at this field's median speed (",
      signif(scale$speed, 2), " m/s) carries a parcel about ", signif(reach, 2),
      " km, and the domain is ", signif(scale$width_km, 2), " by ",
      signif(scale$height_km, 2), " km. ", left_domain, " of ", particles,
      " particles left the velocity field before the window was up."),
    advice = paste0(
      "Shorten integration_days, or fetch a larger bounding box. Either way a ",
      "margin of roughly speed x integration_days is lost along the ", edge,
      " edge, so the usable area is always smaller than the area fetched.")
  )

  env_dat[[name %||% paste0(direction, "_ftle")]] <- result
  env_dat
}

#' Advect particles through a time-varying velocity field
#'
#' Fourth-order Runge-Kutta in lon/lat, with velocities converted from m/s to
#' degrees per day at each evaluation point.
#'
#' @param positions two-column matrix of starting lon/lat
#' @param start_time start time, in days
#' @param duration signed integration length in days; negative runs backwards
#' @param step_days step size in days
#' @param velocity list of two-layer `SpatRaster`s, one per time step
#' @param times numeric time of each raster, in days
#' @return two-column matrix of final lon/lat; `NA` for particles that left the
#'   domain or ran aground
#' @keywords internal
advect <- function(positions, start_time, duration, step_days, velocity, times) {
  n_steps <- max(1, ceiling(abs(duration) / step_days))
  dt <- duration / n_steps
  time <- start_time

  for (step in seq_len(n_steps)) {
    k1 <- velocity_at(positions, time, velocity, times)
    k2 <- velocity_at(positions + 0.5 * dt * k1, time + 0.5 * dt, velocity, times)
    k3 <- velocity_at(positions + 0.5 * dt * k2, time + 0.5 * dt, velocity, times)
    k4 <- velocity_at(positions + dt * k3, time + dt, velocity, times)

    positions <- positions + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    time <- time + dt
  }
  positions
}

#' Velocity in degrees per day at arbitrary positions and time
#'
#' Bilinear in space, linear in time between the two bracketing fields. Sampling
#' the two bracketing rasters at the particle positions and blending the values is
#' equivalent to blending the rasters first, and much cheaper.
#'
#' @param positions two-column matrix of lon/lat
#' @param time time in days
#' @param velocity list of two-layer `SpatRaster`s
#' @param times numeric time of each raster, in days
#' @return two-column matrix of dlon/dt and dlat/dt, in degrees per day
#' @keywords internal
velocity_at <- function(positions, time, velocity, times) {
  blended <- sample_velocity(positions, time, velocity, times)

  # m/s -> degrees/day. Longitude degrees shrink with the cosine of latitude, so
  # the same eastward speed moves a particle through more longitude near the pole.
  seconds_per_day <- 86400
  metres_per_degree <- 111320
  latitude <- positions[, 2]
  cos_lat <- cos(latitude * pi / 180)
  # Guard the pole, where the conversion diverges.
  cos_lat[abs(latitude) > 89.9] <- NA_real_

  cbind(
    blended[, 1] * seconds_per_day / (metres_per_degree * cos_lat),
    blended[, 2] * seconds_per_day / metres_per_degree
  )
}

#' Sample the velocity field at positions and a time
#'
#' @inheritParams velocity_at
#' @return two-column matrix of eastward and northward velocity, in m/s
#' @keywords internal
sample_velocity <- function(positions, time, velocity, times) {
  if (length(times) == 1) {
    return(as.matrix(terra::extract(velocity[[1]], positions, method = "bilinear")))
  }

  # Clamp to the ends of the record rather than extrapolating: a particle that
  # runs past the available fields is held in the last known flow, which is a
  # milder error than inventing velocities.
  time <- min(max(time, min(times)), max(times))
  upper <- findInterval(time, times, all.inside = TRUE) + 1
  lower <- upper - 1

  span <- times[upper] - times[lower]
  weight <- if (span > 0) (time - times[lower]) / span else 0

  before <- as.matrix(terra::extract(velocity[[lower]], positions, method = "bilinear"))
  after <- as.matrix(terra::extract(velocity[[upper]], positions, method = "bilinear"))
  (1 - weight) * before + weight * after
}

#' FTLE from a flow map
#'
#' Differentiates the flow map to get the deformation gradient, forms the
#' Cauchy-Green strain tensor, and returns the exponential separation rate.
#'
#' Positions are converted to metres on an equirectangular plane about the grid
#' centre before differentiating, so the strain tensor is computed in a single
#' consistent metric rather than mixing degrees of longitude and latitude, which
#' have different physical lengths.
#'
#' @param seeds two-column matrix of starting lon/lat, on a regular grid
#' @param final two-column matrix of ending lon/lat
#' @param integration_days integration length in days, unsigned
#' @return numeric vector of FTLE values, in 1/day
#' @keywords internal
ftle_from_flow_map <- function(seeds, final, integration_days) {
  lon <- sort(unique(seeds[, 1]))
  lat <- sort(unique(seeds[, 2]))
  nx <- length(lon)
  ny <- length(lat)

  if (nx < 3 || ny < 3) {
    stop("An FTLE grid needs at least three points in each direction; central ",
         "differences have no interior otherwise.", call. = FALSE)
  }

  metres_per_degree <- 111320
  reference_lat <- mean(lat)
  cos_ref <- cos(reference_lat * pi / 180)

  to_x <- function(longitude) longitude * metres_per_degree * cos_ref
  to_y <- function(latitude) latitude * metres_per_degree

  # Lay the flow map out on the seed lattice so neighbours can be differenced.
  col <- match(seeds[, 1], lon)
  row <- match(seeds[, 2], lat)
  index <- cbind(col, row)

  final_x <- matrix(NA_real_, nx, ny)
  final_y <- matrix(NA_real_, nx, ny)
  final_x[index] <- to_x(final[, 1])
  final_y[index] <- to_y(final[, 2])

  dx <- (lon[2] - lon[1]) * metres_per_degree * cos_ref
  dy <- (lat[2] - lat[1]) * metres_per_degree

  interior_x <- 2:(nx - 1)
  interior_y <- 2:(ny - 1)

  # Deformation gradient by central differences on the seed grid.
  f11 <- (final_x[interior_x + 1, interior_y] - final_x[interior_x - 1, interior_y]) / (2 * dx)
  f12 <- (final_x[interior_x, interior_y + 1] - final_x[interior_x, interior_y - 1]) / (2 * dy)
  f21 <- (final_y[interior_x + 1, interior_y] - final_y[interior_x - 1, interior_y]) / (2 * dx)
  f22 <- (final_y[interior_x, interior_y + 1] - final_y[interior_x, interior_y - 1]) / (2 * dy)

  # Cauchy-Green strain tensor C = F'F, and its largest eigenvalue. For a
  # symmetric 2x2 the eigenvalues are closed-form, so no per-cell eigen() call.
  c11 <- f11^2 + f21^2
  c12 <- f11 * f12 + f21 * f22
  c22 <- f12^2 + f22^2

  trace <- c11 + c22
  gap <- sqrt(pmax((c11 - c22)^2 + 4 * c12^2, 0))
  lambda_max <- (trace + gap) / 2

  # lambda < 1 means net convergence over the window, which is a real result but
  # not an exponential separation rate; clamp so the log stays defined.
  values <- matrix(NA_real_, nx, ny)
  values[interior_x, interior_y] <- log(sqrt(pmax(lambda_max, 1))) / integration_days

  values[index]
}
