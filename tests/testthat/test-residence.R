# Residence time.
#
# Analytic where it can be: with a steady eastward current of known speed, the
# time to cross a box of known width is arithmetic, so the integration can be
# checked against a number rather than against itself.
#
# The case that matters most is censoring. A particle still inside when the
# window ends has a residence of *at least* max_days, and recording it as equal
# to max_days biases every summary downwards. The tests pin both the value and
# the warning, because a silently censored column is the failure that would
# actually reach a model.

steady_field <- function(u, v = 0, months = 1:3,
                         lon = seq(-70, -68, by = 0.25),
                         lat = seq(42, 44, by = 0.25)) {
  frames <- lapply(months, function(m) {
    grid <- expand.grid(x = lon, y = lat)
    grid$YEAR <- 2020L
    grid$MONTH <- as.integer(m)
    grid$DAY <- 1L
    grid$UO <- u
    grid$VO <- v
    grid
  })
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}

# Degrees of longitude covered per day at a given latitude by a current in m/s.
degrees_per_day <- function(speed, latitude) {
  speed * 86400 / (111320 * cos(latitude * pi / 180))
}


test_that("still water never leaves, and says it is censored", {
  env <- steady_field(u = 0, v = 0)
  box <- list(xmin = -69.5, xmax = -68.5, ymin = 42.5, ymax = 43.5)

  expect_warning(result <- residence_time(env, box, max_days = 10), "censored")

  inside <- !is.na(result$forward_residence)
  expect_true(all(result$forward_residence[inside] == 10))
})

test_that("points outside the box get no residence time", {
  env <- steady_field(u = 0.1)
  box <- list(xmin = -69.5, xmax = -69.0, ymin = 42.5, ymax = 43.0)

  suppressWarnings(result <- residence_time(env, box, max_days = 5))

  xy <- sf::st_coordinates(env)
  outside <- xy[, 1] < -69.5 | xy[, 1] > -69.0
  expect_true(all(is.na(result$forward_residence[outside])))
})

test_that("a steady current gives the crossing time the geometry implies", {
  # 0.5 m/s eastward. A particle on the western edge crosses a 1 degree box in
  # the time that speed needs to cover a degree of longitude at this latitude.
  speed <- 0.5
  env <- steady_field(u = speed)
  box <- list(xmin = -69.5, xmax = -68.5, ymin = 42.9, ymax = 43.1)

  step_hours <- 3
  result <- suppressWarnings(
    residence_time(env, box, max_days = 60, step_hours = step_hours)
  )

  xy <- sf::st_coordinates(env)
  west <- abs(xy[, 1] - -69.5) < 1e-9 & abs(xy[, 2] - 43) < 1e-9
  expected <- 1 / degrees_per_day(speed, 43)
  got <- unique(result$forward_residence[west])

  # The box is only checked once a step, so the answer is the first check after
  # the particle actually left: at or above the analytic time, by less than one
  # step.
  expect_gte(got, expected)
  expect_lt(got, expected + step_hours / 24)
})

test_that("a finer step resolves the crossing more closely", {
  speed <- 0.5
  env <- steady_field(u = speed)
  box <- list(xmin = -69.5, xmax = -68.5, ymin = 42.9, ymax = 43.1)
  expected <- 1 / degrees_per_day(speed, 43)

  west_value <- function(step_hours) {
    r <- suppressWarnings(
      residence_time(env, box, max_days = 60, step_hours = step_hours))
    xy <- sf::st_coordinates(env)
    west <- abs(xy[, 1] - -69.5) < 1e-9 & abs(xy[, 2] - 43) < 1e-9
    unique(r$forward_residence[west])
  }

  expect_lt(abs(west_value(1) - expected), abs(west_value(12) - expected))
})

test_that("residence falls the further downstream a particle starts", {
  env <- steady_field(u = 0.3)
  box <- list(xmin = -69.5, xmax = -68.5, ymin = 42.9, ymax = 43.1)

  result <- suppressWarnings(residence_time(env, box, max_days = 60))

  xy <- sf::st_coordinates(env)
  mid <- abs(xy[, 2] - 43) < 1e-9 & !is.na(result$forward_residence)
  by_lon <- tapply(result$forward_residence[mid], xy[mid, 1], function(z) z[1])
  # Ordered west to east, so residence decreases.
  expect_true(all(diff(by_lon) < 0))
})

test_that("a faster current empties the box sooner", {
  box <- list(xmin = -69.5, xmax = -68.5, ymin = 42.9, ymax = 43.1)
  slow <- suppressWarnings(
    residence_time(steady_field(u = 0.1), box, max_days = 60))
  fast <- suppressWarnings(
    residence_time(steady_field(u = 0.5), box, max_days = 60))

  expect_gt(mean(slow$forward_residence, na.rm = TRUE),
            mean(fast$forward_residence, na.rm = TRUE))
})

test_that("backward asks how long water has already been there", {
  # Eastward flow: going backwards, a particle on the eastern edge has been in
  # the box a long time and one on the western edge has just arrived.
  env <- steady_field(u = 0.3)
  box <- list(xmin = -69.5, xmax = -68.5, ymin = 42.9, ymax = 43.1)

  result <- suppressWarnings(
    residence_time(env, box, max_days = 60, direction = "backward"))

  expect_true("backward_residence" %in% names(result))
  xy <- sf::st_coordinates(env)
  mid <- abs(xy[, 2] - 43) < 1e-9 & !is.na(result$backward_residence)
  by_lon <- tapply(result$backward_residence[mid], xy[mid, 1], function(z) z[1])
  # Ordered west to east, so time already spent increases.
  expect_true(all(diff(by_lon) > 0))
})

test_that("censoring is capped at max_days and warned about", {
  # Slow current, short window: nothing gets out.
  env <- steady_field(u = 0.01)
  box <- list(xmin = -69.5, xmax = -68.5, ymin = 42.9, ymax = 43.1)

  expect_warning(result <- residence_time(env, box, max_days = 2), "censored")
  expect_equal(max(result$forward_residence, na.rm = TRUE), 2)
})

test_that("a long enough window is not censored and does not warn", {
  env <- steady_field(u = 0.5)
  box <- list(xmin = -69.5, xmax = -69.0, ymin = 42.9, ymax = 43.1)

  expect_no_warning(result <- residence_time(env, box, max_days = 60))
  expect_lt(max(result$forward_residence, na.rm = TRUE), 60)
})

test_that("an empty box is reported rather than returning all NA quietly", {
  env <- steady_field(u = 0.1)
  box <- list(xmin = 10, xmax = 11, ymin = 10, ymax = 11)

  expect_warning(result <- residence_time(env, box), "nothing was released")
  expect_true(all(is.na(result$forward_residence)))
})

test_that("the default box is the extent of the data", {
  env <- steady_field(u = 0.05)

  default <- suppressWarnings(residence_time(env, max_days = 5))
  explicit <- suppressWarnings(
    residence_time(env, box = extent_box(env), max_days = 5))

  expect_equal(default$forward_residence, explicit$forward_residence)
})

test_that("the column can be named", {
  env <- steady_field(u = 0.1)
  box <- list(xmin = -69.5, xmax = -68.5, ymin = 42.9, ymax = 43.1)

  result <- suppressWarnings(
    residence_time(env, box, max_days = 5, name = "gom_retention"))

  expect_true("gom_retention" %in% names(result))
})

test_that("inside_box excludes NA positions", {
  box <- list(xmin = -70, xmax = -68, ymin = 42, ymax = 44)
  positions <- rbind(c(-69, 43), c(NA, 43), c(-69, NA), c(-60, 43))
  expect_equal(inside_box(positions, box), c(TRUE, FALSE, FALSE, FALSE))
})

test_that("bad arguments are rejected", {
  env <- steady_field(u = 0.1)
  expect_error(residence_time(env, max_days = 0), "positive")
  expect_error(residence_time(env, step_hours = 0), "positive")
  expect_error(residence_time(env, box = list(xmin = 1)), "xmin")
})
