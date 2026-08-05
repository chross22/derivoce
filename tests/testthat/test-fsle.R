test_that("FSLE recovers the analytic exponent of a saddle flow", {
  # A pure strain field u = a*x, v = -a*y separates parcels at exactly rate `a`,
  # so the FSLE must equal `a` regardless of the separation threshold chosen.
  # That invariance is the point of the diagnostic, and this pins it numerically.
  rate_per_day <- 0.1
  metres_per_degree <- 111320
  cos_ref <- cos(42.5 * pi / 180)

  env <- make_flow(function(lon, lat) {
    x <- (lon + 69) * metres_per_degree * cos_ref
    y <- (lat - 42.5) * metres_per_degree
    cbind(rate_per_day * x / 86400, -rate_per_day * y / 86400)
  }, lon = seq(-69.4, -68.6, by = 0.1), lat = seq(42.2, 42.8, by = 0.1))

  result <- fsle(env, final_separation = 20, max_days = 120, step_hours = 6)
  values <- result$backward_fsle[!is.na(result$backward_fsle)]

  expect_gt(length(values), 0)
  expect_equal(median(values), rate_per_day, tolerance = 0.1)
})

test_that("the exponent is insensitive to the separation threshold in a steady saddle", {
  # In a uniform-strain flow the separation rate is scale-independent, so
  # asking for 15 km and 30 km must give the same exponent. In a real flow they
  # would differ - that is the scale selectivity FSLE exists to provide.
  rate_per_day <- 0.1
  metres_per_degree <- 111320
  cos_ref <- cos(42.5 * pi / 180)

  env <- make_flow(function(lon, lat) {
    x <- (lon + 69) * metres_per_degree * cos_ref
    y <- (lat - 42.5) * metres_per_degree
    cbind(rate_per_day * x / 86400, -rate_per_day * y / 86400)
  }, lon = seq(-69.4, -68.6, by = 0.1), lat = seq(42.2, 42.8, by = 0.1))

  near <- fsle(env, final_separation = 15, max_days = 120)
  far <- fsle(env, final_separation = 30, max_days = 120)

  expect_equal(median(near$backward_fsle, na.rm = TRUE),
               median(far$backward_fsle, na.rm = TRUE), tolerance = 0.15)
})

test_that("a stronger strain field separates parcels faster", {
  metres_per_degree <- 111320
  cos_ref <- cos(42.5 * pi / 180)

  saddle <- function(rate) {
    make_flow(function(lon, lat) {
      x <- (lon + 69) * metres_per_degree * cos_ref
      y <- (lat - 42.5) * metres_per_degree
      cbind(rate * x / 86400, -rate * y / 86400)
    }, lon = seq(-69.4, -68.6, by = 0.1), lat = seq(42.2, 42.8, by = 0.1))
  }

  weak <- fsle(saddle(0.05), final_separation = 20, max_days = 120)
  strong <- fsle(saddle(0.15), final_separation = 20, max_days = 120)

  expect_lt(median(weak$backward_fsle, na.rm = TRUE),
            median(strong$backward_fsle, na.rm = TRUE))
})

test_that("parcels that never separate far enough are NA, not slow", {
  # A uniform translation never pulls neighbours apart, so "time to separate by
  # 50 km" has no answer. Returning max_days instead would report a small but
  # non-zero exponent, inventing structure where there is none.
  env <- make_flow(function(lon, lat) cbind(rep(0.2, length(lon)), rep(0, length(lon))))

  result <- fsle(env, final_separation = 50, max_days = 10)

  expect_true(all(is.na(result$backward_fsle)))
})

test_that("direction is recorded in the default column name", {
  metres_per_degree <- 111320
  cos_ref <- cos(42.5 * pi / 180)
  env <- make_flow(function(lon, lat) {
    x <- (lon + 69) * metres_per_degree * cos_ref
    cbind(0.1 * x / 86400, rep(0, length(lon)))
  })

  expect_true("backward_fsle" %in% names(fsle(env, final_separation = 20)))
  expect_true("forward_fsle" %in%
                names(fsle(env, final_separation = 20, direction = "forward")))
  expect_true("scale20" %in%
                names(fsle(env, final_separation = 20, name = "scale20")))
})

test_that("the initial separation defaults to about one grid cell", {
  env <- make_flow(function(lon, lat) cbind(rep(0, length(lon)), rep(0, length(lon))),
                   lon = seq(-70, -68, by = 0.1), lat = seq(42, 43, by = 0.1))
  seeds <- sf::st_coordinates(env[env$MONTH == 1, ])

  # The smaller of the two cell dimensions, so both companion particles start
  # within one cell. At 42.5N, 0.1 degrees of longitude is 8.2 km while 0.1
  # degrees of latitude is 11.1 km, so longitude sets the separation.
  expect_equal(default_separation(seeds), 8.21, tolerance = 0.01)
})

test_that("an initial separation at or beyond the target is rejected", {
  env <- make_flow(function(lon, lat) cbind(rep(0.1, length(lon)), rep(0, length(lon))))

  expect_error(fsle(env, final_separation = 10, initial_separation = 10),
               "must be smaller")
  expect_error(fsle(env, final_separation = 10, initial_separation = 20),
               "must be smaller")
})

test_that("invalid settings are rejected", {
  env <- make_flow(function(lon, lat) cbind(rep(0.1, length(lon)), rep(0, length(lon))))

  expect_error(fsle(env, final_separation = -1), "must all be positive")
  expect_error(fsle(env, final_separation = 20, max_days = 0), "must all be positive")
  expect_error(fsle(env, u = "nope"), "nope")
})
