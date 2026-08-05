make_velocity <- function(u_fun, v_fun = function(lon, lat, month) 0,
                          months = 1:6) {
  lon <- seq(-70, -69, by = 0.25)
  lat <- seq(42, 43, by = 0.25)
  frames <- list()
  for (m in months) {
    grid <- expand.grid(x = lon, y = lat)
    grid$UO <- u_fun(grid$x, grid$y, m)
    grid$VO <- v_fun(grid$x, grid$y, m)
    grid$YEAR <- 2020
    grid$MONTH <- m
    grid$DAY <- 1L
    frames[[length(frames) + 1]] <- grid
  }
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}

test_that("current_speed is the magnitude of the velocity vector", {
  env <- make_velocity(function(lon, lat, month) rep(3, length(lon)),
                       function(lon, lat, month) rep(4, length(lon)))

  result <- current_speed(env)

  expect_equal(result$speed, rep(5, nrow(env)))
})

test_that("a steady current has zero eddy kinetic energy", {
  # Strong but unvarying: high total kinetic energy, no eddy energy. A version
  # that forgot to subtract the mean would report 0.5 * 0.5^2 here.
  env <- make_velocity(function(lon, lat, month) rep(0.5, length(lon)))

  result <- eke(env)

  expect_true(all(abs(result$EKE) < 1e-12))
})

test_that("EKE matches the analytic value for a known fluctuation", {
  # Alternates +0.2 and -0.2 about a mean of zero, so every anomaly is 0.2 and
  # EKE is 0.5 * 0.2^2 = 0.02 everywhere.
  env <- make_velocity(function(lon, lat, month) {
    rep(ifelse(month %% 2 == 0, 0.2, -0.2), length(lon))
  })

  result <- eke(env)

  expect_equal(result$EKE, rep(0.02, nrow(env)))
})

test_that("EKE uses both velocity components", {
  env <- make_velocity(
    function(lon, lat, month) rep(ifelse(month %% 2 == 0, 0.3, -0.3), length(lon)),
    function(lon, lat, month) rep(ifelse(month %% 2 == 0, 0.4, -0.4), length(lon))
  )

  result <- eke(env)

  # 0.5 * (0.3^2 + 0.4^2)
  expect_equal(result$EKE, rep(0.125, nrow(env)))
})

test_that("the mean is removed per location, not globally", {
  # Half the domain flows east, half flows west, both steadily. A single
  # domain-wide mean would leave a large spurious anomaly everywhere; a
  # per-location mean correctly gives zero.
  env <- make_velocity(function(lon, lat, month) ifelse(lon > -69.5, 0.5, -0.5))

  result <- eke(env)

  expect_true(all(abs(result$EKE) < 1e-12))
})

test_that("climatology removes a repeating seasonal cycle that record does not", {
  # A pure seasonal cycle, identical every year. Against the record mean it
  # registers as energetic; against a per-month climatology it vanishes, which
  # is the entire distinction between the two references.
  env <- make_velocity(function(lon, lat, month) rep(sin(month / 12 * 2 * pi), length(lon)),
                       months = 1:12)
  env2 <- env
  env2$YEAR <- 2021
  both <- rbind(env, env2)

  against_record <- eke(both)
  against_climatology <- eke(both, reference = "climatology")

  expect_gt(mean(against_record$EKE), 0.1)
  expect_true(all(abs(against_climatology$EKE) < 1e-12))
})

test_that("a rolling reference keeps only variability faster than its window", {
  # A slow linear drift with no fast fluctuation. A short rolling mean tracks
  # the drift and removes it; the record mean cannot, and reports the drift as
  # eddy energy.
  env <- make_velocity(function(lon, lat, month) rep(month * 0.1, length(lon)),
                       months = 1:12)

  against_record <- eke(env)
  against_rolling <- eke(env, reference = 3)

  expect_gt(mean(against_record$EKE), mean(against_rolling$EKE))
})

test_that("invalid references are rejected", {
  env <- make_velocity(function(lon, lat, month) rep(0.1, length(lon)))

  expect_error(eke(env, reference = 1), "at least 2")
  expect_error(eke(env, reference = "seasonal"),
               "must be \"record\", \"climatology\", or a positive integer")
})

test_that("missing velocity columns are reported", {
  env <- make_velocity(function(lon, lat, month) rep(0.1, length(lon)))

  expect_error(eke(env, u = "nope"), "nope")
  expect_error(current_speed(env, v = "nope"), "Available")
})
