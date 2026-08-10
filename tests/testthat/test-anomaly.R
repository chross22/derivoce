# Cell-wise anomalies.
#
# The failure this guards hardest against is the silent one: a monthly
# climatology built from a single year returns a column of exact zeroes, which
# looks like a perfectly ordinary covariate and carries no information at all.

anomaly_field <- function(values, years, months,
                          lon = c(-70, -69.5), lat = c(43, 43.5)) {
  grid <- expand.grid(x = lon, y = lat)
  frames <- list()
  for (y in years) {
    for (m in months) {
      f <- grid
      f$YEAR <- y
      f$MONTH <- m
      f$DAY <- 1L
      f$SST <- values(f$x, f$y, y, m)
      frames[[length(frames) + 1]] <- f
    }
  }
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}


test_that("a cell that never changes has no anomaly", {
  env <- anomaly_field(function(lon, lat, y, m) lon * 0 + 10,
                       years = 2000:2004, months = 1:3)

  # Constant in time, which the package warns about for a temporal operation.
  expect_warning(result <- cell_anomaly(env, "SST"), "[Ss]tatic")

  expect_true(all(abs(result$SST_anom) < 1e-12))
})

test_that("the anomaly is the departure from that cell's own mean", {
  # One cell warms by 1 degree a year; the climatology is the middle year.
  env <- anomaly_field(function(lon, lat, y, m) (y - 2002) + 10,
                       years = 2000:2004, months = 1L)

  result <- cell_anomaly(env, "SST")

  expect_equal(sort(unique(round(result$SST_anom, 10))), c(-2, -1, 0, 1, 2))
})

test_that("each cell is compared against itself, not against its neighbours", {
  # A strong west-east background gradient with identical variability on top.
  env <- anomaly_field(function(lon, lat, y, m) 100 * lon + (y - 2002),
                       years = 2000:2004, months = 1L)

  result <- cell_anomaly(env, "SST")

  # The background is removed entirely: only the shared departure survives.
  by_year <- tapply(result$SST_anom, env$YEAR, function(z) unique(round(z, 10)))
  expect_equal(as.numeric(unlist(by_year)), c(-2, -1, 0, 1, 2))
})

test_that("climatology compares a month against the same month", {
  # A large seasonal cycle and one anomalous January.
  env <- anomaly_field(function(lon, lat, y, m) {
    seasonal <- 10 * m
    seasonal + ifelse(y == 2003 & m == 1, 5, 0)
  }, years = 2000:2004, months = 1:3)

  result <- cell_anomaly(env, "SST", reference = "climatology")

  odd <- env$YEAR == 2003 & env$MONTH == 1
  expect_true(all(result$SST_anom[odd] > 3))
  expect_true(all(abs(result$SST_anom[!odd]) < 2))
})

test_that("record leaves the seasonal cycle in and climatology removes it", {
  env <- anomaly_field(function(lon, lat, y, m) 10 * m,
                       years = 2000:2004, months = 1:3)

  clim <- cell_anomaly(env, "SST", reference = "climatology")
  rec <- cell_anomaly(env, "SST", reference = "record")

  expect_true(all(abs(clim$SST_anom) < 1e-12))
  expect_gt(stats::sd(rec$SST_anom), 5)
})

test_that("standardising gives a z-score", {
  env <- anomaly_field(function(lon, lat, y, m) (y - 2002) + 10,
                       years = 2000:2004, months = 1L)

  result <- cell_anomaly(env, "SST", standardize = TRUE)

  # Values -2..2 have sd sqrt(2.5); the z-scores are those divided by it.
  expect_equal(sort(unique(round(result$SST_z, 8))),
               round(sort(c(-2, -1, 0, 1, 2) / stats::sd(-2:2)), 8))
  expect_equal(mean(result$SST_z), 0)
})

test_that("standardising makes places of different variability comparable", {
  # West varies by 1 degree, east by 10. The raw anomalies differ tenfold; the
  # z-scores do not.
  env <- anomaly_field(function(lon, lat, y, m) {
    amplitude <- ifelse(lon < -69.7, 1, 10)
    amplitude * (y - 2002)
  }, years = 2000:2004, months = 1L)

  raw <- cell_anomaly(env, "SST")
  z <- cell_anomaly(env, "SST", standardize = TRUE)
  west <- sf::st_coordinates(env)[, 1] < -69.7

  expect_gt(stats::sd(raw$SST_anom[!west]) / stats::sd(raw$SST_anom[west]), 9)
  expect_equal(stats::sd(z$SST_z[!west]), stats::sd(z$SST_z[west]))
})

test_that("a single year of monthly data warns instead of returning zeroes", {
  env <- anomaly_field(function(lon, lat, y, m) 10 * m + lon,
                       years = 2000L, months = 1:12)

  expect_warning(result <- cell_anomaly(env, "SST"), "exactly zero")
  # The warning is the point: the numbers really are all zero.
  expect_true(all(abs(result$SST_anom) < 1e-12))
})

test_that("the thin-climatology warning names the way out", {
  env <- anomaly_field(function(lon, lat, y, m) 10 * m + lon,
                       years = 2000L, months = 1:12)

  expect_warning(cell_anomaly(env, "SST"), "record")
})

test_that("record over a single year is fine and does not warn", {
  env <- anomaly_field(function(lon, lat, y, m) 10 * m + lon,
                       years = 2000L, months = 1:12)

  expect_no_warning(result <- cell_anomaly(env, "SST", reference = "record"))
  expect_gt(stats::sd(result$SST_anom), 1)
})

test_that("missing values do not poison a cell's mean", {
  env <- anomaly_field(function(lon, lat, y, m) (y - 2002) + 10,
                       years = 2000:2004, months = 1L)
  env$SST[env$YEAR == 2004] <- NA

  result <- cell_anomaly(env, "SST")

  # Mean of -2..1 is -0.5, so the 2000 value sits 1.5 below it.
  expect_equal(unique(round(result$SST_anom[env$YEAR == 2000], 10)), -1.5)
  expect_true(all(is.na(result$SST_anom[env$YEAR == 2004])))
})

test_that("a suffix can be chosen", {
  env <- anomaly_field(function(lon, lat, y, m) (y - 2002),
                       years = 2000:2004, months = 1L)

  result <- cell_anomaly(env, "SST", suffix = "_departure")

  expect_true("SST_departure" %in% names(result))
})

test_that("several covariates are handled in one call", {
  env <- anomaly_field(function(lon, lat, y, m) (y - 2002),
                       years = 2000:2004, months = 1L)
  env$SSS <- 32 + (env$YEAR - 2002) * 0.1

  result <- cell_anomaly(env, c("SST", "SSS"))

  expect_true(all(c("SST_anom", "SSS_anom") %in% names(result)))
  expect_equal(mean(result$SSS_anom), 0)
})

test_that("a missing covariate is an error", {
  env <- anomaly_field(function(lon, lat, y, m) (y - 2002),
                       years = 2000:2004, months = 1L)
  expect_error(cell_anomaly(env, "NOPE"), "NOPE")
})


test_that("detrending removes a linear trend the plain anomaly leaves in", {
  # Ten years warming steadily, with no seasonal cycle.
  env <- anomaly_field(function(lon, lat, y, m) 0.5 * (y - 2000),
                       years = 2000:2009, months = 1L)

  plain <- cell_anomaly(env, "SST", reference = "record")
  flat <- cell_anomaly(env, "SST", reference = "record", detrend = TRUE)

  # The plain anomaly is the trend, and correlates almost perfectly with year.
  expect_gt(abs(stats::cor(plain$SST_anom, env$YEAR)), 0.99)
  # Detrended, what is left is a ten-thousandth of what was there. Not machine
  # zero, because the series is linear in *year* while the trend is fitted
  # against elapsed days, and leap years make those spacings uneven.
  expect_lt(stats::sd(flat$SST_anom), stats::sd(plain$SST_anom) / 1000)
})

test_that("a warm year stands out once the trend is gone", {
  # A steady warming plus one genuinely anomalous year. Against the raw
  # anomaly the last years all look warm; against the detrended one, only the
  # anomalous year does.
  env <- anomaly_field(function(lon, lat, y, m) {
    0.5 * (y - 2000) + ifelse(y == 2003, 4, 0)
  }, years = 2000:2009, months = 1L)

  flat <- cell_anomaly(env, "SST", reference = "record", detrend = TRUE)

  odd <- env$YEAR == 2003
  expect_equal(which.max(flat$SST_anom), which(odd)[1])
  expect_gt(min(flat$SST_anom[odd]), max(flat$SST_anom[!odd]))
})

test_that("detrending with a climatology drops the seasonal cycle as well", {
  env <- anomaly_field(function(lon, lat, y, m) {
    0.3 * (y - 2000) + 5 * sin(2 * pi * m / 12)
  }, years = 2000:2009, months = 1:12)

  both <- cell_anomaly(env, "SST", reference = "climatology", detrend = TRUE)
  kept <- cell_anomaly(env, "SST", reference = "record", detrend = TRUE)

  # Trend and cycle both gone.
  expect_lt(stats::sd(both$SST_anom), 0.05)
  # Trend gone, cycle still there.
  expect_gt(stats::sd(kept$SST_anom), 3)
})

test_that("standardising still works when detrending", {
  env <- anomaly_field(function(lon, lat, y, m) {
    0.5 * (y - 2000) + ifelse(y %% 2 == 0, 1, -1)
  }, years = 2000:2009, months = 1L)

  z <- cell_anomaly(env, "SST", reference = "record", detrend = TRUE,
                    standardize = TRUE)

  expect_true("SST_z" %in% names(z))
  expect_equal(mean(z$SST_z), 0, tolerance = 1e-8)
  expect_equal(stats::sd(z$SST_z), 1, tolerance = 0.2)
})

test_that("detrending is off by default, so behaviour is unchanged", {
  env <- anomaly_field(function(lon, lat, y, m) 0.5 * (y - 2000),
                       years = 2000:2009, months = 1L)

  expect_equal(cell_anomaly(env, "SST", reference = "record")$SST_anom,
               cell_anomaly(env, "SST", reference = "record",
                            detrend = FALSE)$SST_anom)
})
