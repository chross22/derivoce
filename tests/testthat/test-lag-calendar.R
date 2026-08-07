# Lags by calendar unit rather than by position in the series. The value at each
# step encodes its own date, so the source a lag drew from is identifiable.

stamped <- function(years, months, days = 1L) {
  frames <- list()
  for (y in years) for (m in months) for (d in days) {
    grid <- expand.grid(x = seq(-70, -69, by = 0.5), y = seq(42, 43, by = 0.5))
    grid$SST <- y * 10000 + m * 100 + d
    grid$YEAR <- y; grid$MONTH <- m; grid$DAY <- as.integer(d)
    frames[[length(frames) + 1]] <- grid
  }
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}

stamp <- function(y, m, d = 1) y * 10000 + m * 100 + d

at <- function(env, column, y, m, d = 1) {
  unique(env[[column]][env$YEAR == y & env$MONTH == m & env$DAY == d])
}

test_that("by = 'step' is unchanged, and remains the default", {
  env <- stamped(2020, 1:6)

  result <- lag_covariate(env, "SST")

  expect_true("SST_lag1" %in% names(result))
  expect_equal(at(result, "SST_lag1", 2020, 5), stamp(2020, 4))
  expect_true(is.na(at(result, "SST_lag1", 2020, 1)))
})

test_that("month lags count calendar months and cross the year boundary", {
  env <- stamped(2020:2021, 1:12)

  result <- lag_covariate(env, "SST", n = 3, by = "month")

  expect_true("SST_lag3month" %in% names(result))
  expect_equal(at(result, "SST_lag3month", 2020, 5), stamp(2020, 2))
  expect_equal(at(result, "SST_lag3month", 2021, 2), stamp(2020, 11))
  expect_true(is.na(at(result, "SST_lag3month", 2020, 2)))
})

test_that("year lags hold the calendar month fixed", {
  env <- stamped(2020:2021, 1:12)

  result <- lag_covariate(env, "SST", n = 1, by = "year")

  expect_equal(at(result, "SST_lag1year", 2021, 7), stamp(2020, 7))
  expect_true(is.na(at(result, "SST_lag1year", 2020, 7)))
})

test_that("day lags use real dates", {
  env <- stamped(2020, 1, days = 1:10)

  result <- lag_covariate(env, "SST", n = 7, by = "day")

  expect_equal(at(result, "SST_lag7day", 2020, 1, 9), stamp(2020, 1, 2))
  expect_true(is.na(at(result, "SST_lag7day", 2020, 1, 3)))
})

test_that("a gap in the record separates 'step' from 'month'", {
  # This is the whole reason the argument exists. April is missing, so May's
  # previous *step* is March: a one-step lag is quietly a two-month one, and
  # nothing in the output says so.
  env <- stamped(2020, c(1, 2, 3, 5, 6))

  by_step <- lag_covariate(env, "SST", n = 1, by = "step")
  by_month <- lag_covariate(env, "SST", n = 1, by = "month")

  expect_equal(at(by_step, "SST_lag1", 2020, 5), stamp(2020, 3))
  expect_true(is.na(at(by_month, "SST_lag1month", 2020, 5)))
  # Where there is no gap the two agree.
  expect_equal(at(by_month, "SST_lag1month", 2020, 3), stamp(2020, 2))
})

test_that("a vector of lags gives an autoregressive set", {
  env <- stamped(2018:2022, 1:12)

  result <- lag_covariate(env, "SST", n = 1:3, by = "year")

  expect_true(all(c("SST_lag1year", "SST_lag2year", "SST_lag3year") %in%
                    names(result)))
  expect_equal(at(result, "SST_lag1year", 2022, 7), stamp(2021, 7))
  expect_equal(at(result, "SST_lag2year", 2022, 7), stamp(2020, 7))
  expect_equal(at(result, "SST_lag3year", 2022, 7), stamp(2019, 7))
  # The deepest lag runs out first.
  expect_true(is.na(at(result, "SST_lag3year", 2019, 7)))
})

test_that("suffixes are named per lag, or supplied per lag", {
  env <- stamped(2020:2022, 1:12)

  expect_true(all(c("SST_a", "SST_b") %in%
                    names(lag_covariate(env, "SST", n = 1:2, by = "year",
                                        suffix = c("_a", "_b")))))
  expect_error(lag_covariate(env, "SST", n = 1:3, by = "year", suffix = "_one"),
               "one entry per lag")
})

test_that("n is validated", {
  env <- stamped(2020, 1:6)

  expect_error(lag_covariate(env, "SST", n = 0), "at least 1")
  expect_error(lag_covariate(env, "SST", n = 1.5), "whole number")
  expect_error(lag_covariate(env, "SST", n = NA), "whole number")
})

test_that("the source-step map is computable from the step table alone", {
  # The calendar arithmetic is separable from the data, which is what makes it
  # testable without building a grid for every case.
  steps <- data.frame(YEAR = c(2020, 2020, 2020, 2020),
                      MONTH = c(1, 2, 3, 5), DAY = 1L)

  expect_equal(lag_source_step(steps, 1, "step"), c(NA, 1, 2, 3))
  # April is absent, so May (row 4) has no one-month predecessor.
  expect_equal(lag_source_step(steps, 1, "month"), c(NA, 1, 2, NA))
  # Two months before May is March, which is row 3 rather than row 2. Position
  # and calendar have already come apart by here.
  expect_equal(lag_source_step(steps, 2, "month"), c(NA, NA, 1, 3))
  expect_equal(lag_source_step(steps, 1, "year"), rep(NA_integer_, 4))
})
