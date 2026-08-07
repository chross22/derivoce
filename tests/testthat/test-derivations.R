# Analytic checks on the non-Lagrangian derivations: fields whose answer is known
# in closed form, so a wrong result is caught by value rather than by looking
# plausible. Complements the per-function test files, which cover shape and
# argument handling.

# ---- eddy kinetic energy ---------------------------------------------------

test_that("EKE against the record mean recovers the analytic value", {
  # u alternates +a, -a between two steps, v is zero. The per-location mean is
  # therefore exactly 0, the anomaly is +/-a, and EKE = (a^2 + 0)/2.
  a <- 0.3
  env <- make_flow(function(lon, lat) cbind(rep(0, length(lon)), rep(0, length(lon))),
                   months = 1:2)
  env$UO <- ifelse(env$MONTH == 1, a, -a)
  env$VO <- 0

  result <- eke(env)

  expect_equal(unique(result$EKE), a^2 / 2)
})

test_that("EKE against a climatology removes a repeating seasonal cycle", {
  # The same two-month cycle in two consecutive years. Every value equals its own
  # calendar-month mean, so the anomaly is identically zero. This is the property
  # that distinguishes "climatology" from "record", and it must hold exactly.
  a <- 0.3
  one_year <- make_flow(function(lon, lat) cbind(rep(0, length(lon)), rep(0, length(lon))),
                        months = 1:2)
  one_year$UO <- ifelse(one_year$MONTH == 1, a, -a)
  one_year$VO <- 0
  env <- sf::st_as_sf(rbind(transform(one_year, YEAR = 2010),
                            transform(one_year, YEAR = 2011)))

  expect_true(all(abs(eke(env, reference = "climatology")$EKE) < 1e-12))
  # The same data against the record mean must NOT be zero, or the test above
  # would pass for a function that always returned zero.
  expect_true(all(eke(env, reference = "record")$EKE > 0))
})

test_that("current speed is the Pythagorean magnitude", {
  env <- make_flow(function(lon, lat) cbind(rep(3, length(lon)), rep(4, length(lon))),
                   months = 1)

  expect_true(all(abs(current_speed(env)$speed - 5) < 1e-12))
})

# ---- temporal derivations --------------------------------------------------

test_that("per = 'day' divides by real calendar days, not by step", {
  # January to February is 31 days; February to March is 29 in 2020. A per-step
  # rate cannot tell them apart, which is the whole reason `per` exists.
  env <- make_env(function(lon, lat, year, month) month * 10, months = 1:3)

  per_step <- temporal_gradient(env, "SST", per = "step")$SST_tgrad
  per_day <- temporal_gradient(env, "SST", per = "day")$SST_tgrad

  expect_true(all(abs(stats::na.omit(per_step) - 10) < 1e-9))

  jan_feb <- as.numeric(as.Date("2020-02-01") - as.Date("2020-01-01"))
  feb_mar <- as.numeric(as.Date("2020-03-01") - as.Date("2020-02-01"))
  expect_equal(sort(unique(round(stats::na.omit(per_day), 9))),
               sort(round(c(10 / jan_feb, 10 / feb_mar), 9)))
})

test_that("the integration window resets at the year boundary, or does not", {
  # The covariate equals the month number, so each window's total is a small sum
  # that can be written down by hand. A covariate of 1 everywhere would be
  # simpler still, but it is static in time and would (correctly) trip the
  # degeneracy warning, which belongs to a different test.
  env <- make_env(function(lon, lat, year, month) month, years = 2020:2021, months = 1:3)
  at <- function(x, y, m) x[env$YEAR == y & env$MONTH == m][1]

  yearly <- integrate_covariate(env, "SST", window = "year")$SST_int
  expect_equal(at(yearly, 2020, 3), 1 + 2 + 3)
  expect_equal(at(yearly, 2021, 1), 1)   # reset at the year boundary

  cumulative <- integrate_covariate(env, "SST", window = "all")$SST_int
  expect_equal(at(cumulative, 2021, 1), 1 + 2 + 3 + 1)   # no reset

  rolling <- integrate_covariate(env, "SST", window = 2)$SST_int
  expect_equal(at(rolling, 2021, 3), 2 + 3)   # trailing two steps
})

test_that("lags survive a shuffled time step", {
  # Matching by row order would silently give every row some other cell's value.
  env <- make_env(function(lon, lat, year, month) lon + lat + month, months = 1:3)
  shuffled <- env
  rows <- which(shuffled$MONTH == 2)
  set.seed(1)
  shuffled[rows, ] <- shuffled[rows, ][sample(length(rows)), ]

  plain <- lag_covariate(env, "SST")
  mixed <- lag_covariate(shuffled, "SST")

  stamp <- function(z) paste(location_key(z), z$YEAR, z$MONTH, z$DAY)
  expect_equal(plain$SST_lag1, mixed$SST_lag1[match(stamp(plain), stamp(mixed))])
  expect_gt(sum(is.finite(plain$SST_lag1)), 0)
})

# ---- vertical gradient -----------------------------------------------------

test_that("vertical gradient is a difference, or a per-metre rate with depth", {
  env <- make_env(function(lon, lat, year, month) 12)
  env$BOTT <- 4
  env$DEPTH <- 200

  total <- vertical_gradient(env)
  expect_true(all(total$SST_BOTT_vgrad == 8))

  rate <- vertical_gradient(env, depth = "DEPTH")
  expect_true(all(abs(rate$SST_BOTT_vgrad - 8 / 200) < 1e-12))
})

test_that("zero depth gives NA rather than Inf", {
  env <- make_env(function(lon, lat, year, month) 12)
  env$BOTT <- 4
  env$DEPTH <- 0

  expect_true(all(is.na(vertical_gradient(env, depth = "DEPTH")$SST_BOTT_vgrad)))
})

# ---- distances -------------------------------------------------------------

test_that("distance to a contour matches the analytic distance, to one cell", {
  # A linear depth ramp puts the 100 m contour at exactly lon = -68, so the
  # distance from any point is a known function of longitude.
  #
  # The tolerance is one grid cell, and that is a property of the method rather
  # than slack: terra::distance() measures to the centre of the nearest cell
  # flagged as on-contour, and the flagged cells straddle the true crossing. So
  # the error is bounded by the cell size in absolute terms, which makes the
  # RELATIVE error large close in and small far away.
  lon <- seq(-70, -68, by = 0.05)
  lat <- seq(42, 43, by = 0.05)
  grid <- expand.grid(x = lon, y = lat)
  grid$DEPTH <- 50 + 50 * (grid$x + 70) / 2
  grid$YEAR <- 2020; grid$MONTH <- 1; grid$DAY <- 1L
  env <- sf::st_as_sf(grid, coords = c("x", "y"), crs = 4326)

  result <- distance_to_contour(env, "DEPTH", levels = 100)
  xy <- sf::st_coordinates(result)
  expected <- abs(xy[, 1] + 68) * 111.32 * cos(xy[, 2] * pi / 180)

  cell_km <- 0.05 * 111.32 * cos(42.5 * pi / 180)
  got <- result$DEPTH_dist_100
  expect_true(max(abs(got - expected), na.rm = TRUE) < 1.1 * cell_km)
  expect_lt(min(got, na.rm = TRUE), 1)
})

test_that("distance units scale exactly, and the isobath wrapper agrees", {
  lon <- seq(-70, -68, by = 0.1)
  lat <- seq(42, 43, by = 0.1)
  grid <- expand.grid(x = lon, y = lat)
  grid$DEPTH <- 50 + 50 * (grid$x + 70) / 2
  grid$YEAR <- 2020; grid$MONTH <- 1; grid$DAY <- 1L
  env <- sf::st_as_sf(grid, coords = c("x", "y"), crs = 4326)

  km <- distance_to_contour(env, "DEPTH", levels = 100)$DEPTH_dist_100
  m <- distance_to_contour(env, "DEPTH", levels = 100, per = "m")$DEPTH_dist_100
  expect_equal(m, km * 1000)

  isobath <- distance_to_isobath(env, levels = 100)$isobath_dist_100
  expect_equal(isobath, km)
})

test_that("front distance is zero on a front and grows away from it", {
  lon <- seq(-70, -68, by = 0.05)
  lat <- seq(42, 43, by = 0.05)
  grid <- expand.grid(x = lon, y = lat)
  grid$SST <- ifelse(grid$x > -69, 15, 10)   # a step, so the front is at -69
  grid$YEAR <- 2020; grid$MONTH <- 1; grid$DAY <- 1L
  env <- sf::st_as_sf(grid, coords = c("x", "y"), crs = 4326)

  result <- distance_to_front(env, "SST", quantile = 0.95)
  xy <- sf::st_coordinates(result)
  distance <- result$SST_front_dist

  on_front <- abs(xy[, 1] + 69) < 0.06
  far_away <- abs(xy[, 1] + 69) > 0.9
  expect_lt(stats::median(distance[on_front], na.rm = TRUE), 10)
  expect_gt(stats::median(distance[far_away], na.rm = TRUE), 50)
})

# ---- float32 coordinates ----------------------------------------------------

test_that("a grid with float32 coordinates is accepted, as Copernicus supplies", {
  # Copernicus stores longitude and latitude as float32. Near 67 degrees that
  # resolves to about 8e-6, so a nominally uniform 1/12-degree grid arrives with
  # spacings varying by roughly 1e-4 of a cell. Every synthetic fixture in this
  # suite uses exact seq() doubles and cannot show this, which is how a tolerance
  # that rejected every real download survived until someone ran one.
  round_trip <- function(x) as.numeric(as.vector(writeBin(
    as.numeric(x), raw(), size = 4) |> readBin("numeric", n = length(x), size = 4)))

  lon <- round_trip(seq(-67.17, -64.83, by = 1 / 12))
  lat <- round_trip(seq(42, 44.17, by = 1 / 12))

  # The quantisation is real and larger than the old 1e-6 tolerance allowed.
  expect_gt(max(abs(diff(lon) - stats::median(diff(lon)))) /
              stats::median(diff(lon)), 1e-6)
  expect_true(is_regular(lon))
  expect_true(is_regular(lat))

  grid <- expand.grid(x = lon, y = lat)
  grid$SST <- grid$x + grid$y
  grid$YEAR <- 2010; grid$MONTH <- 1; grid$DAY <- 1L
  env <- sf::st_as_sf(grid, coords = c("x", "y"), crs = 4326)

  # terra::rast(type = "xyz") applies its own regularity check, so the raster is
  # built from the lattice dimensions instead. This is what proves the two agree.
  expect_no_error(result <- horizontal_gradient(env, "SST"))
  expect_gt(sum(is.finite(result$SST_grad)), 0)
})

test_that("genuinely irregular spacing is still rejected", {
  # The looser tolerance must not let real irregularity through. Both of these
  # are off by orders of magnitude more than float32 noise.
  expect_false(is_regular(c(0, 1, 2, 4, 5)))          # a missing column
  expect_false(is_regular(sort(runif(20))))            # scattered
  expect_true(is_regular(seq(0, 1, by = 0.1)))         # exact
})

test_that("values land in the right cells after the explicit rasterization", {
  # Building the raster by index rather than from coordinates could transpose or
  # flip it, and a gradient would still look plausible. A field that varies only
  # with latitude pins the orientation.
  lon <- seq(-70, -68, by = 0.25)
  lat <- seq(42, 43, by = 0.25)
  grid <- expand.grid(x = lon, y = lat)
  grid$SST <- grid$y * 10
  grid$YEAR <- 2020; grid$MONTH <- 1; grid$DAY <- 1L
  env <- sf::st_as_sf(grid, coords = c("x", "y"), crs = 4326)

  result <- horizontal_gradient(env, "SST", components = TRUE)
  interior_rows <- interior(env, margin = 0.25)

  # All signal northward, none eastward, and the value is the analytic one.
  expect_true(all(abs(result$SST_grad_x[interior_rows]) < 1e-9))
  expect_true(all(abs(result$SST_grad_y[interior_rows] - 10 / 111.32) < 1e-6))
})
