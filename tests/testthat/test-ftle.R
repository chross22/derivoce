test_that("a motionless field has zero FTLE", {
  # Nothing moves, so nothing separates and the strain tensor is the identity.
  env <- make_flow(function(lon, lat) cbind(rep(0, length(lon)), rep(0, length(lon))))

  result <- ftle(env, integration_days = 5, step_hours = 12)
  values <- result$backward_ftle[!is.na(result$backward_ftle)]

  expect_gt(length(values), 0)
  expect_true(all(values < 1e-10))
})

test_that("a constant m/s flow gives a small but real FTLE from meridian convergence", {
  # A constant eastward speed in m/s is NOT a uniform flow in lon/lat: meridians
  # converge toward the pole, so a parcel at 43N sweeps through more longitude
  # per second than one at 42N. That differential is genuine shear, and the
  # resulting small non-zero FTLE is the correct answer, not numerical noise.
  # Across 42-43N the spread is cos(42)/cos(43) - 1, about 1.6%.
  env <- make_flow(function(lon, lat) cbind(rep(0.1, length(lon)), rep(0, length(lon))))

  result <- ftle(env, integration_days = 5, step_hours = 12)
  values <- result$backward_ftle[!is.na(result$backward_ftle)]

  expect_true(all(values > 0))
  # Small: orders of magnitude below the strain rates of real ocean fronts,
  # which run 0.01-0.1 per day.
  expect_true(all(values < 1e-3))
})

test_that("a saddle flow recovers its analytic exponent", {
  # A pure strain field u = a*x, v = -a*y separates neighbouring particles at
  # exactly rate `a` along one axis, so the FTLE must equal `a` regardless of
  # integration time. This is the sharpest available check: it pins the
  # numerical value, not just its sign or ordering.
  rate_per_day <- 0.05
  metres_per_degree <- 111320
  centre_lon <- -69
  centre_lat <- 42.5
  cos_ref <- cos(centre_lat * pi / 180)

  # Convert the strain rate from per-day to the m/s field it implies.
  env <- make_flow(function(lon, lat) {
    x <- (lon - centre_lon) * metres_per_degree * cos_ref
    y <- (lat - centre_lat) * metres_per_degree
    cbind(rate_per_day * x / 86400, -rate_per_day * y / 86400)
  })

  result <- ftle(env, integration_days = 10, step_hours = 6)
  values <- result$backward_ftle[!is.na(result$backward_ftle)]

  expect_equal(median(values), rate_per_day, tolerance = 0.02)
})

test_that("forward and backward integration both detect a symmetric saddle", {
  # A steady saddle is time-reversible, so its attracting and repelling
  # structures are mirror images and the two directions give the same
  # magnitude. On an unsteady field they would differ.
  rate_per_day <- 0.05
  metres_per_degree <- 111320
  cos_ref <- cos(42.5 * pi / 180)

  env <- make_flow(function(lon, lat) {
    x <- (lon + 69) * metres_per_degree * cos_ref
    y <- (lat - 42.5) * metres_per_degree
    cbind(rate_per_day * x / 86400, -rate_per_day * y / 86400)
  })

  backward <- ftle(env, integration_days = 10, direction = "backward")
  forward <- ftle(env, integration_days = 10, direction = "forward")

  expect_true("backward_ftle" %in% names(backward))
  expect_true("forward_ftle" %in% names(forward))
  expect_equal(median(backward$backward_ftle, na.rm = TRUE),
               median(forward$forward_ftle, na.rm = TRUE),
               tolerance = 0.05)
})

test_that("direction is recorded in the default column name", {
  env <- make_flow(function(lon, lat) cbind(rep(0.05, length(lon)), rep(0, length(lon))))

  expect_true("backward_ftle" %in% names(ftle(env, integration_days = 2)))
  expect_true("forward_ftle" %in%
                names(ftle(env, integration_days = 2, direction = "forward")))
  expect_true("lcs" %in% names(ftle(env, integration_days = 2, name = "lcs")))
})

test_that("a stronger strain field gives a larger exponent", {
  metres_per_degree <- 111320
  cos_ref <- cos(42.5 * pi / 180)

  saddle <- function(rate) {
    make_flow(function(lon, lat) {
      x <- (lon + 69) * metres_per_degree * cos_ref
      y <- (lat - 42.5) * metres_per_degree
      cbind(rate * x / 86400, -rate * y / 86400)
    })
  }

  weak <- ftle(saddle(0.02), integration_days = 10)
  strong <- ftle(saddle(0.08), integration_days = 10)

  expect_lt(median(weak$backward_ftle, na.rm = TRUE),
            median(strong$backward_ftle, na.rm = TRUE))
})

test_that("grid boundary cells are NA", {
  env <- make_flow(function(lon, lat) cbind(rep(0.1, length(lon)), rep(0, length(lon))))

  result <- ftle(env, integration_days = 2)
  coords <- sf::st_coordinates(result)
  on_edge <- coords[, 1] == min(coords[, 1]) | coords[, 1] == max(coords[, 1]) |
    coords[, 2] == min(coords[, 2]) | coords[, 2] == max(coords[, 2])

  # Central differences have no interior at the edge, so those cells cannot
  # have a value.
  expect_true(all(is.na(result$backward_ftle[on_edge])))
  expect_false(all(is.na(result$backward_ftle[!on_edge])))
})

test_that("invalid integration settings are rejected", {
  env <- make_flow(function(lon, lat) cbind(rep(0.1, length(lon)), rep(0, length(lon))))

  # Direction is an argument, so a negative duration is a mistake rather than a
  # second way of saying "backward".
  expect_error(ftle(env, integration_days = -5), "must be positive")
  expect_error(ftle(env, integration_days = 5, step_hours = 0), "must be positive")
  expect_error(ftle(env, integration_days = 1, step_hours = 48),
               "no longer than the integration window")
})

test_that("missing velocity columns are reported", {
  env <- make_flow(function(lon, lat) cbind(rep(0.1, length(lon)), rep(0, length(lon))))

  expect_error(ftle(env, u = "missing_u"), "missing_u")
  expect_error(ftle(env, u = "missing_u"), "Available")
})

test_that("a grid too small to differentiate is rejected", {
  env <- make_flow(function(lon, lat) cbind(rep(0.1, length(lon)), rep(0, length(lon))),
                   lon = c(-70, -69.9), lat = c(42, 42.1))

  expect_error(ftle(env, integration_days = 2), "at least three points")
})
