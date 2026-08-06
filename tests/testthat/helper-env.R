# A regular lon/lat grid across several time steps, in the shape
# datamatch::accessEnvDat() returns: one sf POINT row per (grid point, time step)
# with YEAR/MONTH/DAY columns.
#
# `fun(lon, lat, year, month)` supplies the covariate value, so a test can plant
# a field whose derivative is known analytically and check the computed one
# against it.
make_env <- function(fun = function(lon, lat, year, month) lon + lat,
                     lon = seq(-70, -66, by = 0.5),
                     lat = seq(41, 44, by = 0.5),
                     years = 2020,
                     months = 1,
                     var = "SST") {
  frames <- list()
  for (year in years) {
    for (month in months) {
      grid <- expand.grid(x = lon, y = lat)
      grid[[var]] <- fun(grid$x, grid$y, year, month)
      grid$YEAR <- year
      grid$MONTH <- month
      grid$DAY <- 1L
      frames[[length(frames) + 1]] <- grid
    }
  }
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}

# Interior points only. Central differences are undefined on the boundary, where
# one of the two neighbours is missing, so edge cells are legitimately NA.
interior <- function(env, margin = 0.5) {
  coords <- sf::st_coordinates(env)
  which(coords[, 1] > min(coords[, 1]) + margin - 1e-9 &
          coords[, 1] < max(coords[, 1]) - margin + 1e-9 &
          coords[, 2] > min(coords[, 2]) + margin - 1e-9 &
          coords[, 2] < max(coords[, 2]) - margin + 1e-9)
}

# TRUE if evaluating `expr` raised a warning. Lets a test assert that a check
# stayed quiet, which is as much a part of its contract as when it speaks.
warned <- function(expr) {
  tryCatch({ expr; FALSE }, warning = function(w) TRUE)
}

# Velocity fields as sf points, in the shape accessEnvDat() returns. `fun`
# returns a two-column matrix of (u, v) in m/s for the given lon/lat.
make_flow <- function(fun,
                      lon = seq(-70, -68, by = 0.05),
                      lat = seq(42, 43, by = 0.05),
                      months = 1:2) {
  frames <- list()
  for (m in months) {
    grid <- expand.grid(x = lon, y = lat)
    uv <- fun(grid$x, grid$y)
    grid$UO <- uv[, 1]
    grid$VO <- uv[, 2]
    grid$YEAR <- 2020
    grid$MONTH <- m
    grid$DAY <- 1L
    frames[[length(frames) + 1]] <- grid
  }
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}
