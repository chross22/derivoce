test_that("a uniform field has zero gradient", {
  env <- make_env(function(lon, lat, year, month) rep(10, length(lon)))

  result <- horizontal_gradient(env, "SST")

  expect_true(all(abs(result$SST_grad[interior(env)]) < 1e-9))
})

test_that("gradient magnitude matches an analytic north-south ramp", {
  # 1 degree C per degree of latitude. A degree of latitude is 111.32 km
  # everywhere, so the answer is 1/111.32 degrees C per km at every point,
  # independent of longitude and latitude.
  env <- make_env(function(lon, lat, year, month) lat)

  result <- horizontal_gradient(env, "SST", components = TRUE)
  rows <- interior(env)

  expect_equal(result$SST_grad[rows], rep(1 / 111.32, length(rows)), tolerance = 1e-4)
  # The ramp is purely northward, so the eastward component must vanish.
  expect_true(all(abs(result$SST_grad_x[rows]) < 1e-9))
  expect_equal(result$SST_grad_y[rows], rep(1 / 111.32, length(rows)), tolerance = 1e-4)
})

test_that("an east-west ramp is scaled by the cosine of latitude", {
  # 1 degree C per degree of longitude. Unlike latitude, the distance in a
  # degree of longitude shrinks toward the pole, so the same temperature change
  # is spread over less distance and the gradient is steeper further north.
  # Ignoring that - the easy mistake - would give a constant here.
  env <- make_env(function(lon, lat, year, month) lon)

  result <- horizontal_gradient(env, "SST")
  coords <- sf::st_coordinates(env)
  rows <- interior(env)

  expected <- 1 / (111.32 * cos(coords[rows, 2] * pi / 180))
  expect_equal(result$SST_grad[rows], expected, tolerance = 1e-3)
  # Concretely: the gradient at 44 N is measurably steeper than at 41 N.
  north <- result$SST_grad[rows][which.max(coords[rows, 2])]
  south <- result$SST_grad[rows][which.min(coords[rows, 2])]
  expect_gt(north, south)
})

test_that("the distance unit changes the magnitude by a factor of 1000", {
  env <- make_env(function(lon, lat, year, month) lat)

  per_km <- horizontal_gradient(env, "SST")$SST_grad
  per_m <- horizontal_gradient(env, "SST", per = "m")$SST_grad
  rows <- interior(env)

  expect_equal(per_km[rows], per_m[rows] * 1000, tolerance = 1e-9)
})

test_that("gradients are computed per time step, not across the whole record", {
  # The field is flat within each month but jumps between months. A gradient
  # that pooled time steps would see that jump as spatial structure.
  env <- make_env(function(lon, lat, year, month) rep(month * 10, length(lon)),
                  months = 1:3)

  result <- horizontal_gradient(env, "SST")

  expect_true(all(abs(result$SST_grad[interior(env)]) < 1e-9))
  expect_equal(nrow(result), nrow(env))
})

test_that("scattered points are rejected rather than silently interpolated", {
  env <- make_env()
  # Jitter the coordinates so they no longer form a regular lattice.
  scattered <- sf::st_as_sf(
    data.frame(x = runif(50, -70, -66), y = runif(50, 41, 44), SST = runif(50),
               YEAR = 2020, MONTH = 1, DAY = 1),
    coords = c("x", "y"), crs = 4326
  )

  expect_error(horizontal_gradient(scattered, "SST"), "not on a regular")
})

test_that("an unknown covariate lists the available ones", {
  env <- make_env()

  expect_error(horizontal_gradient(env, "NOPE"), "NOPE")
  expect_error(horizontal_gradient(env, "NOPE"), "Available")
})

test_that("vertical gradient is the surface-minus-bottom difference", {
  env <- make_env()
  env$BOTT <- env$SST - 4

  result <- vertical_gradient(env, "SST", "BOTT")

  expect_equal(result$SST_BOTT_vgrad, rep(4, nrow(env)))
})

test_that("vertical gradient can be expressed per metre of depth", {
  env <- make_env()
  env$BOTT <- env$SST - 4
  env$DEPTH <- 200

  result <- vertical_gradient(env, "SST", "BOTT", depth = "DEPTH")

  expect_equal(result$SST_BOTT_vgrad, rep(4 / 200, nrow(env)))
})

test_that("zero or negative depth yields NA rather than Inf", {
  env <- make_env()
  env$BOTT <- env$SST - 4
  env$DEPTH <- c(0, rep(200, nrow(env) - 1))

  result <- vertical_gradient(env, "SST", "BOTT", depth = "DEPTH")

  expect_true(is.na(result$SST_BOTT_vgrad[1]))
  expect_false(any(is.infinite(result$SST_BOTT_vgrad)))
})

test_that("temporal gradient measures change between consecutive steps", {
  # Warms by 2 degrees a month, everywhere.
  env <- make_env(function(lon, lat, year, month) rep(month * 2, length(lon)),
                  months = 1:4)

  result <- temporal_gradient(env, "SST")

  # The first step has no predecessor.
  expect_true(all(is.na(result$SST_tgrad[result$MONTH == 1])))
  expect_equal(result$SST_tgrad[result$MONTH == 2], rep(2, sum(result$MONTH == 2)))
  expect_equal(result$SST_tgrad[result$MONTH == 4], rep(2, sum(result$MONTH == 4)))
})

test_that("temporal gradient can be expressed per day", {
  env <- make_env(function(lon, lat, year, month) rep(month * 31, length(lon)),
                  months = 1:2)

  result <- temporal_gradient(env, "SST", per = "day")

  # Jan 1 to Feb 1 is 31 days, and the field rose by 31, so 1 per day.
  expect_equal(result$SST_tgrad[result$MONTH == 2], rep(1, sum(result$MONTH == 2)))
})
