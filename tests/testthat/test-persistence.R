# How often a place is frontal.
#
# The measure exists to separate a persistent feature from a procession of
# transient ones, so the central test is a domain containing both: a front that
# never moves and a front that is somewhere different every step. An
# implementation that averaged distance to the nearest front, or that counted
# any step with a front anywhere, would rank the two the same.

# A step whose SST has a sharp jump at a given longitude.
front_at <- function(edge) {
  function(lon, lat) ifelse(lon < edge, 4, 12)
}

front_field <- function(edges, lon = seq(-70, -68, by = 0.1),
                        lat = seq(43, 43.4, by = 0.1)) {
  frames <- lapply(seq_along(edges), function(i) {
    grid <- expand.grid(x = lon, y = lat)
    grid$YEAR <- 2020L
    grid$MONTH <- as.integer(i)
    grid$DAY <- 1L
    grid$SST <- front_at(edges[i])(grid$x, grid$y)
    grid
  })
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}

longitudes <- function(env) round(sf::st_coordinates(env)[, 1], 6)


test_that("a stationary front scores higher than a wandering one", {
  # West: the same edge every step. East: a different edge every step.
  lon <- seq(-70, -68, by = 0.1)
  edges_west <- rep(-69.5, 6)
  frames <- lapply(1:6, function(i) {
    grid <- expand.grid(x = lon, y = seq(43, 43.4, by = 0.1))
    grid$YEAR <- 2020L
    grid$MONTH <- as.integer(i)
    grid$DAY <- 1L
    # A fixed step at -69.5 and a moving one in the eastern half.
    moving <- -68.9 + 0.1 * i
    grid$SST <- ifelse(grid$x < -69.5, 4, 12) +
      ifelse(grid$x < moving, 0, 6)
    grid
  })
  env <- sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)

  result <- front_frequency(env, "SST", scope = "step")

  x <- longitudes(env)
  freq <- result$SST_front_freq
  stationary <- mean(freq[abs(x - -69.5) < 0.06], na.rm = TRUE)
  wandering <- mean(freq[x > -68.9], na.rm = TRUE)
  expect_gt(stationary, wandering)
})

test_that("a front that never moves is frontal in every step", {
  env <- front_field(rep(-69, 5))

  result <- front_frequency(env, "SST", scope = "step")

  x <- longitudes(env)
  at_front <- abs(x - -69) < 0.06
  expect_equal(max(result$SST_front_freq[at_front], na.rm = TRUE), 1)
})

test_that("frequencies are proportions", {
  env <- front_field(c(-69, -69, -69, -68.5, -68.5))

  result <- front_frequency(env, "SST", scope = "step")
  values <- result$SST_front_freq[!is.na(result$SST_front_freq)]

  expect_true(all(values >= 0 & values <= 1))
})

test_that("a place frontal in three of five steps scores three fifths", {
  # The front sits at -69 for three steps and at -68.5 for two.
  env <- front_field(c(-69, -69, -69, -68.5, -68.5))

  result <- front_frequency(env, "SST", scope = "step", quantile = 0.9)

  x <- longitudes(env)
  at_first <- abs(x - -69) < 0.06
  expect_equal(max(result$SST_front_freq[at_first], na.rm = TRUE), 0.6,
               tolerance = 1e-8)
})

test_that("the whole record gives one static value per cell", {
  env <- front_field(c(-69, -69, -68.5))

  result <- front_frequency(env, "SST", scope = "step")

  # The same cell has the same value in every step.
  key <- paste(longitudes(env), round(sf::st_coordinates(env)[, 2], 6))
  per_cell <- tapply(result$SST_front_freq, key,
                     function(z) length(unique(round(z[!is.na(z)], 10))))
  expect_true(all(per_cell <= 1))
})

test_that("a window makes it vary through the record", {
  env <- front_field(c(-69, -69, -69, -68.5, -68.5, -68.5))

  result <- front_frequency(env, "SST", scope = "step", n = 2)

  x <- longitudes(env)
  at_first <- abs(x - -69) < 0.06
  early <- result$SST_front_freq[at_first & result$MONTH == 3]
  late <- result$SST_front_freq[at_first & result$MONTH == 6]

  # Frontal in both of the last two steps early on, in neither by the end.
  expect_equal(max(early, na.rm = TRUE), 1)
  expect_equal(max(late, na.rm = TRUE), 0)
})

test_that("an undefined gradient leaves the denominator", {
  # The outermost ring has no central difference, so it is never counted as
  # frontal and never counted as not frontal either.
  env <- front_field(rep(-69, 4))

  result <- front_frequency(env, "SST", scope = "step")

  x <- longitudes(env)
  edge <- abs(x - min(x)) < 1e-9
  expect_true(all(is.na(result$SST_front_freq[edge])))
})

test_that("an explicit threshold is honoured", {
  env <- front_field(rep(-69, 4))

  # Far above any gradient present, so nothing is frontal.
  none <- front_frequency(env, "SST", threshold = 1e6)
  expect_equal(max(none$SST_front_freq, na.rm = TRUE), 0)

  # Far below, so everywhere the field actually changes is frontal. Flat water
  # has a gradient of exactly zero and remains not frontal at any positive
  # threshold, so this does not reach every cell.
  loose <- front_frequency(env, "SST", threshold = 1e-9)
  expect_equal(max(loose$SST_front_freq, na.rm = TRUE), 1)
  expect_equal(min(loose$SST_front_freq, na.rm = TRUE), 0)

  # A step field has only two gradient levels, zero and the jump, so every
  # threshold between them selects the same cells. Thresholds separate cells on
  # a field with structure, which is checked on the quantile path elsewhere.
  default <- front_frequency(env, "SST")
  expect_gte(mean(loose$SST_front_freq > 0, na.rm = TRUE),
             mean(default$SST_front_freq > 0, na.rm = TRUE))
})

test_that("the column can be named", {
  env <- front_field(rep(-69, 3))
  result <- front_frequency(env, "SST", name = "thermal_persistence")
  expect_true("thermal_persistence" %in% names(result))
})

test_that("bad arguments are rejected", {
  env <- front_field(rep(-69, 3))
  expect_error(front_frequency(env, "SST", threshold = -1), "positive")
  expect_error(front_frequency(env, "SST", quantile = 1), "between 0 and 1")
  expect_error(front_frequency(env, "SST", n = 0), "at least 1")
  expect_error(front_frequency(env, "NOPE"), "NOPE")
})
