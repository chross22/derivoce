# Pulling a broadcast per-step index back out as a series.
#
# The refusal is the interesting part. Collapsing a column that varies across
# the grid would return one arbitrary cell's value and look entirely plausible,
# so a map has to be rejected rather than quietly averaged.

series_field <- function(n_steps = 4, lon = c(-70, -69.5), lat = c(43, 43.5)) {
  grid <- expand.grid(x = lon, y = lat)
  frames <- lapply(seq_len(n_steps), function(i) {
    f <- grid
    f$YEAR <- 2020L
    f$MONTH <- as.integer(i)
    f$DAY <- 1L
    f$index <- i * 10          # constant within a step
    f$map <- i * 10 + f$x      # varies across the grid
    f
  })
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}


test_that("an index collapses to one row per time step", {
  env <- series_field(4)

  out <- index_series(env, "index")

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 4)
  expect_equal(out$index, c(10, 20, 30, 40))
})

test_that("the result carries the date columns and is not an sf object", {
  env <- series_field(3)

  out <- index_series(env, "index")

  expect_true(all(c("YEAR", "MONTH", "DAY") %in% names(out)))
  expect_false(inherits(out, "sf"))
  expect_false("geometry" %in% names(out))
})

test_that("rows come back in date order", {
  env <- series_field(3)
  env <- env[rev(seq_len(nrow(env))), ]

  out <- index_series(env, "index")

  expect_equal(out$MONTH, 1:3)
  expect_equal(out$index, c(10, 20, 30))
})

test_that("a column that varies within a step is refused", {
  env <- series_field(3)

  expect_error(index_series(env, "map"), "vary within a time step")
  # And the message says what to do instead.
  expect_error(index_series(env, "map"), "box_anomaly")
})

test_that("auto-detection finds the indices and leaves the maps alone", {
  env <- series_field(3)

  out <- index_series(env)

  expect_true("index" %in% names(out))
  expect_false("map" %in% names(out))
})

test_that("several indices come back in one call", {
  env <- series_field(3)
  env$other <- env$MONTH * -2

  out <- index_series(env, c("index", "other"))

  expect_equal(out$index, c(10, 20, 30))
  expect_equal(out$other, c(-2, -4, -6))
})

test_that("an object with no index at all is an error that explains", {
  env <- series_field(3)
  env$index <- NULL

  expect_error(index_series(env), "no index to extract")
})

test_that("a missing column is an error naming what is available", {
  env <- series_field(3)
  expect_error(index_series(env, "NOPE"), "NOPE")
})

test_that("it round-trips a real index computed by the package", {
  env <- series_field(6)
  env$SSS <- 32 + env$MONTH * 0.1 + sf::st_coordinates(env)[, 1] * 0.01
  env <- box_anomaly(env, "SSS",
                     box = list(xmin = -71, xmax = -69, ymin = 42, ymax = 44),
                     reference = "record")

  out <- index_series(env, "SSS_box_anom")

  expect_equal(nrow(out), 6)
  # A record-referenced anomaly is centred, so the series sums to zero.
  expect_equal(sum(out$SSS_box_anom), 0, tolerance = 1e-8)
})

test_that("a step where the index is missing comes back NA, not dropped", {
  env <- series_field(4)
  env$index[env$MONTH == 2] <- NA

  out <- index_series(env, "index")

  expect_equal(nrow(out), 4)
  expect_true(is.na(out$index[out$MONTH == 2]))
  expect_equal(out$index[out$MONTH == 3], 30)
})
