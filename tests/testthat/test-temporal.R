test_that("lag_covariate returns the previous step's value at the same location", {
  env <- make_env(function(lon, lat, year, month) rep(month * 10, length(lon)),
                  months = 1:4)

  result <- lag_covariate(env, "SST")

  expect_true("SST_lag1" %in% names(result))
  expect_true(all(is.na(result$SST_lag1[result$MONTH == 1])))
  expect_equal(result$SST_lag1[result$MONTH == 2], rep(10, sum(result$MONTH == 2)))
  expect_equal(result$SST_lag1[result$MONTH == 4], rep(30, sum(result$MONTH == 4)))
})

test_that("lags of more than one step look further back", {
  env <- make_env(function(lon, lat, year, month) rep(month * 10, length(lon)),
                  months = 1:4)

  result <- lag_covariate(env, "SST", n = 2)

  expect_true("SST_lag2" %in% names(result))
  expect_true(all(is.na(result$SST_lag2[result$MONTH %in% 1:2])))
  expect_equal(result$SST_lag2[result$MONTH == 4], rep(20, sum(result$MONTH == 4)))
})

test_that("lagged values follow location, not row order", {
  # A spatially varying field, so a lag that matched by position rather than by
  # coordinate would return the wrong cell's value.
  env <- make_env(function(lon, lat, year, month) lon * 100 + lat + month,
                  months = 1:2)

  # Shuffle the second month's rows.
  second <- which(env$MONTH == 2)
  env[second, ] <- env[sample(second), ]

  result <- lag_covariate(env, "SST")
  rows <- which(result$MONTH == 2)

  # Each February value is exactly 1 more than its own January value.
  expect_equal(result$SST[rows] - result$SST_lag1[rows], rep(1, length(rows)))
})

test_that("lag_covariate rejects a lag below one", {
  expect_error(lag_covariate(make_env(), "SST", n = 0), "at least 1")
})

test_that("integrate_covariate accumulates from the start of each year", {
  # 1 in January, 2 in February, 3 in March. The running total is 1, 3, 6, and
  # resets when the next year starts.
  env <- make_env(function(lon, lat, year, month) rep(month, length(lon)),
                  years = c(2020, 2021), months = 1:3)

  result <- integrate_covariate(env, "SST")

  expect_equal(result$SST_int[result$YEAR == 2020 & result$MONTH == 1][1], 1)
  expect_equal(result$SST_int[result$YEAR == 2020 & result$MONTH == 3][1], 6)
  # The reset is the point of window = "year": 2021 starts over rather than
  # continuing from 6.
  expect_equal(result$SST_int[result$YEAR == 2021 & result$MONTH == 1][1], 1)
  expect_equal(result$SST_int[result$YEAR == 2021 & result$MONTH == 3][1], 6)
})

test_that("a numeric window gives a rolling total that does not reset", {
  env <- make_env(function(lon, lat, year, month) rep(month, length(lon)),
                  years = c(2020, 2021), months = 1:3)

  result <- integrate_covariate(env, "SST", window = 2)

  # Trailing two steps: Feb 2020 is 1 + 2.
  expect_equal(result$SST_int[result$YEAR == 2020 & result$MONTH == 2][1], 3)
  # Jan 2021 sums across the year boundary: Mar 2020 (3) + Jan 2021 (1).
  expect_equal(result$SST_int[result$YEAR == 2021 & result$MONTH == 1][1], 4)
})

test_that("window = all accumulates over the whole record", {
  env <- make_env(function(lon, lat, year, month) rep(month, length(lon)),
                  years = c(2020, 2021), months = 1:3)

  result <- integrate_covariate(env, "SST", window = "all")

  # 1+2+3 in 2020, then 1+2+3 again in 2021, never resetting.
  expect_equal(result$SST_int[result$YEAR == 2021 & result$MONTH == 3][1], 12)
})

test_that("integrate_covariate rejects an unusable window", {
  env <- make_env()

  expect_error(integrate_covariate(env, "SST", window = 0), "at least 1")
  expect_error(integrate_covariate(env, "SST", window = "decade"),
               "must be \"year\", \"all\", or a positive integer")
})

test_that("a location absent from every contributing step is NA, not zero", {
  env <- make_env(function(lon, lat, year, month) rep(month, length(lon)),
                  months = 1:2)
  env$SST[env$MONTH == 1] <- NA
  env$SST[env$MONTH == 2] <- NA

  result <- integrate_covariate(env, "SST")

  # Summing nothing must not look like a genuine total of zero.
  expect_true(all(is.na(result$SST_int)))
})
