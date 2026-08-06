#' Finite-Size Lyapunov Exponent
#'
#' Measures how *quickly* neighbouring water parcels reach a chosen separation,
#' rather than how far apart they get in a chosen time. Where [ftle()] fixes the
#' clock and measures distance, FSLE fixes the distance and measures the clock:
#'
#' \deqn{\lambda = \frac{\ln(\delta_f / \delta_0)}{\tau}}
#'
#' with \eqn{\delta_0} the starting separation, \eqn{\delta_f} the target, and
#' \eqn{\tau} the time taken to reach it. Units are 1/day.
#'
#' @section Why choose this over FTLE:
#' The difference is **scale selectivity**, and it matters most in a domain whose
#' flow speed varies a lot from place to place.
#'
#' An FTLE map with a fixed integration time resolves fine structure where the
#' flow is fast and only coarse structure where it is slow. Across a shelf with an
#' energetic break current and a sluggish interior, ridge intensity then partly
#' encodes background current speed rather than frontal activity, and a model
#' cannot tell the two apart.
#'
#' FSLE asks the same question everywhere — "how long to separate by
#' \eqn{\delta_f}?" — so results are comparable between energetic and quiet
#' regions, and \eqn{\delta_f} can be set to a scale that means something
#' biologically: a patch size, a predator's search radius, a survey's resolution.
#'
#' The cost is that FSLE has no natural place to encode a biologically meaningful
#' *timescale*. When the relevant duration is known — a retention time, a cohort's
#' accumulation window — [ftle()] with that integration time is the better tool.
#'
#' @section How it is computed:
#' Each grid point is seeded with two companion particles, offset east and north
#' by `initial_separation`. All three are advected together, and the separation of
#' each pair is checked after every step. \eqn{\tau} is the first time either pair
#' reaches `final_separation`.
#'
#' Parcels that never separate that far within `max_days` return `NA`: the
#' question "how long to reach this separation" simply has no answer for them, and
#' substituting `max_days` would report a slow separation rate where there was
#' none at all.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`, on a
#'   regular lon/lat grid, containing eastward and northward velocity columns
#' @param u name of the eastward velocity column, in m/s
#' @param v name of the northward velocity column, in m/s
#' @param final_separation target separation \eqn{\delta_f}, in km. This is the
#'   scale-selectivity knob: it sets the size of structure resolved.
#' @param initial_separation starting separation \eqn{\delta_0}, in km; defaults
#'   to roughly one grid cell
#' @param max_days give up after this long. Must exceed the time a typical parcel
#'   pair needs, or most of the field returns `NA`.
#' @param direction `"backward"` (attracting structures, the default) or
#'   `"forward"` (repelling structures)
#' @param step_hours Runge-Kutta step size, in hours
#' @param name name for the new column
#' @return `env_dat` with an FSLE column added, in 1/day. Parcels that never reach
#'   `final_separation` within `max_days` are `NA`.
#' @examples
#' \dontrun{
#' # Structures at the 50 km scale, comparable across the whole shelf
#' env <- datamatch::accessEnvDat(vars = c("UO", "VO"), ...)
#' env <- fsle(env, final_separation = 50)       # UO and VO are the defaults
#'
#' # A scale matched to a predator's search radius
#' env <- fsle(env, final_separation = 10, max_days = 60)
#' }
#' @references
#' d'Ovidio, F., Fernandez, V., Hernandez-Garcia, E., & Lopez, C. (2004). Mixing
#' structures in the Mediterranean Sea from finite-size Lyapunov exponents.
#' *Geophysical Research Letters*, 31(17).
#' @seealso [ftle()], and `docs/methods.md` for when to prefer which
#' @export
fsle <- function(env_dat, u = "UO", v = "VO", final_separation = 50,
                 initial_separation = NULL, max_days = 60,
                 direction = c("backward", "forward"), step_hours = 6,
                 name = NULL) {
  direction <- match.arg(direction)
  resolve_vars(env_dat, c(u, v))

  if (final_separation <= 0 || max_days <= 0 || step_hours <= 0) {
    stop("final_separation, max_days, and step_hours must all be positive.",
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

  for (i in seq_len(nrow(steps))) {
    rows <- step_rows(env_dat, steps[i, ])
    seeds <- sf::st_coordinates(env_dat[rows, ])

    # Resolved per time step, since the default follows the grid rather than
    # being a fixed distance.
    delta_0 <- initial_separation %||% default_separation(seeds)
    if (delta_0 >= final_separation) {
      stop("initial_separation (", signif(delta_0, 3), " km) must be smaller ",
           "than final_separation (", final_separation, " km); parcels cannot ",
           "separate to a distance they already start beyond.", call. = FALSE)
    }

    elapsed <- separation_time(seeds, delta_0, final_separation,
                               start_time = step_times[i], sign = sign,
                               max_days = max_days,
                               step_days = step_hours / 24,
                               velocity = velocity, times = step_times)

    result[rows] <- log(final_separation / delta_0) / elapsed
  }

  env_dat[[name %||% paste0(direction, "_fsle")]] <- result
  env_dat
}

#' Default initial separation for a seeded grid
#'
#' One grid cell, in km, measured at the grid's mean latitude.
#'
#' @param seeds two-column matrix of lon/lat
#' @return separation in km
#' @keywords internal
default_separation <- function(seeds) {
  lon <- sort(unique(seeds[, 1]))
  lat <- sort(unique(seeds[, 2]))
  reference_lat <- mean(lat)
  metres_per_degree <- 111320

  dx <- (lon[2] - lon[1]) * metres_per_degree * cos(reference_lat * pi / 180)
  dy <- (lat[2] - lat[1]) * metres_per_degree
  min(dx, dy) / 1000
}

#' Time for a parcel pair to reach a target separation
#'
#' Advects each seed together with an eastward and a northward companion,
#' checking after every step whether either pair has separated far enough.
#'
#' @param seeds two-column matrix of lon/lat
#' @param delta_0 initial separation, in km
#' @param delta_f target separation, in km
#' @param start_time start time in days
#' @param sign `-1` to integrate backwards, `1` forwards
#' @param max_days give-up time
#' @param step_days step size in days
#' @param velocity list of two-layer `SpatRaster`s
#' @param times numeric time of each raster, in days
#' @return numeric vector of times in days; `NA` where the target was never met
#' @keywords internal
separation_time <- function(seeds, delta_0, delta_f, start_time, sign, max_days,
                            step_days, velocity, times) {
  metres_per_degree <- 111320
  reference_lat <- mean(seeds[, 2])
  offset_lon <- (delta_0 * 1000) / (metres_per_degree * cos(reference_lat * pi / 180))
  offset_lat <- (delta_0 * 1000) / metres_per_degree

  centre <- seeds
  east <- cbind(seeds[, 1] + offset_lon, seeds[, 2])
  north <- cbind(seeds[, 1], seeds[, 2] + offset_lat)

  elapsed <- rep(NA_real_, nrow(seeds))
  n_steps <- ceiling(max_days / step_days)
  dt <- sign * step_days
  time <- start_time

  # Only parcels still short of the target are advected. Without this, every
  # parcel is integrated for the full max_days even after its answer is known,
  # which dominates the cost: in a real field most pairs separate early, and the
  # long tail is a small minority.
  active <- seq_len(nrow(seeds))

  for (step in seq_len(n_steps)) {
    centre[active, ] <- rk4_step(centre[active, , drop = FALSE], time, dt, velocity, times)
    east[active, ] <- rk4_step(east[active, , drop = FALSE], time, dt, velocity, times)
    north[active, ] <- rk4_step(north[active, , drop = FALSE], time, dt, velocity, times)
    time <- time + dt

    separation <- pmax(
      pair_distance(centre[active, , drop = FALSE], east[active, , drop = FALSE], reference_lat),
      pair_distance(centre[active, , drop = FALSE], north[active, , drop = FALSE], reference_lat)
    )

    # A parcel that has drifted off the grid has no velocity and cannot be
    # followed further, so it is retired unanswered rather than carried as NA.
    lost <- is.na(separation)
    reached <- !lost & separation >= delta_f
    elapsed[active[reached]] <- step * step_days

    active <- active[!reached & !lost]
    if (length(active) == 0) break
  }
  elapsed
}

#' One fourth-order Runge-Kutta step
#'
#' @param positions two-column matrix of lon/lat
#' @param time current time in days
#' @param dt signed step size in days
#' @param velocity list of two-layer `SpatRaster`s
#' @param times numeric time of each raster, in days
#' @return two-column matrix of updated lon/lat
#' @keywords internal
rk4_step <- function(positions, time, dt, velocity, times) {
  k1 <- velocity_at(positions, time, velocity, times)
  k2 <- velocity_at(positions + 0.5 * dt * k1, time + 0.5 * dt, velocity, times)
  k3 <- velocity_at(positions + 0.5 * dt * k2, time + 0.5 * dt, velocity, times)
  k4 <- velocity_at(positions + dt * k3, time + dt, velocity, times)

  positions + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
}

#' Distance between paired positions, in km
#'
#' Equirectangular about a reference latitude, which is accurate for the small
#' separations FSLE deals with and consistent with the metric used elsewhere here.
#'
#' @param a,b two-column matrices of lon/lat
#' @param reference_lat latitude to evaluate the longitude scaling at
#' @return numeric vector of distances in km
#' @keywords internal
pair_distance <- function(a, b, reference_lat) {
  metres_per_degree <- 111320
  dx <- (b[, 1] - a[, 1]) * metres_per_degree * cos(reference_lat * pi / 180)
  dy <- (b[, 2] - a[, 2]) * metres_per_degree
  sqrt(dx^2 + dy^2) / 1000
}
