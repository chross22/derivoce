# Rolling summaries over a trailing window.
#
# The step-versus-calendar distinction is the substance here, and it only shows
# up on a record with a gap in it: until then the two agree exactly, which is
# what makes getting it wrong so easy to miss.

rolling_field <- function(values, years = 2020L, months = seq_along(values),
                          lon = c(-70, -69.5), lat = 43) {
  grid <- expand.grid(x = lon, y = lat)
  years <- rep_len(years, length(values))
  frames <- lapply(seq_along(values), function(i) {
    f <- grid
    f$YEAR <- as.integer(years[i])
    f$MONTH <- as.integer(months[i])
    f$DAY <- 1L
    f$SST <- values[i]
    f
  })
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}

by_month <- function(result, column) {
  vapply(split(result[[column]], result$MONTH), function(z) unique(z)[1],
         numeric(1))
}


test_that("the window is trailing and includes the current step", {
  env <- rolling_field(c(1, 2, 3, 4, 5))

  result <- rolling_covariate(env, "SST", n = 3, stat = "mean")

  # March averages January, February and March.
  expect_equal(unname(by_month(result, "SST_mean3")[3]), 2)
  expect_equal(unname(by_month(result, "SST_mean3")[5]), 4)
})

test_that("early steps summarise the short window they have", {
  env <- rolling_field(c(1, 2, 3, 4, 5))

  result <- rolling_covariate(env, "SST", n = 3, stat = "mean")

  expect_equal(unname(by_month(result, "SST_mean3")[1]), 1)
  expect_equal(unname(by_month(result, "SST_mean3")[2]), 1.5)
})

test_that("min_obs withholds a summary of too short a window", {
  env <- rolling_field(c(1, 2, 3, 4, 5))

  result <- rolling_covariate(env, "SST", n = 3, stat = "mean", min_obs = 3)

  expect_true(all(is.na(result$SST_mean3[result$MONTH %in% 1:2])))
  expect_equal(unname(by_month(result, "SST_mean3")[3]), 2)
})

test_that("each statistic computes what it says", {
  env <- rolling_field(c(1, 2, 3, 10, 5))

  result <- rolling_covariate(env, "SST", n = 3,
                              stat = c("mean", "min", "max", "sum", "median",
                                       "range"))

  at4 <- function(col) unname(by_month(result, col)[4])
  expect_equal(at4("SST_mean3"), mean(c(2, 3, 10)))
  expect_equal(at4("SST_min3"), 2)
  expect_equal(at4("SST_max3"), 10)
  expect_equal(at4("SST_sum3"), 15)
  expect_equal(at4("SST_median3"), 3)
  expect_equal(at4("SST_range3"), 8)
})

test_that("a standard deviation needs two values", {
  env <- rolling_field(c(1, 2, 3, 4, 5))

  result <- rolling_covariate(env, "SST", n = 3, stat = "sd")

  expect_true(is.na(unname(by_month(result, "SST_sd3")[1])))
  expect_equal(unname(by_month(result, "SST_sd3")[3]), stats::sd(1:3))
})

test_that("a window of one returns the value itself", {
  env <- rolling_field(c(4, 7, 1))

  result <- rolling_covariate(env, "SST", n = 1, stat = "mean")

  expect_equal(result$SST_mean1, result$SST)
})

test_that("step and calendar windows agree on a complete record", {
  env <- rolling_field(c(1, 2, 3, 4, 5))

  by_step <- rolling_covariate(env, "SST", n = 3, by = "step", stat = "mean")
  by_cal <- rolling_covariate(env, "SST", n = 3, by = "month", stat = "mean")

  expect_equal(by_step$SST_mean3, by_cal$SST_mean3month)
})

test_that("they disagree once the record has a gap, in the documented way", {
  # Monthly series missing April.
  env <- rolling_field(c(1, 2, 3, 5, 6), months = c(1, 2, 3, 5, 6))

  by_step <- rolling_covariate(env, "SST", n = 3, by = "step", stat = "mean")
  by_cal <- rolling_covariate(env, "SST", n = 3, by = "month", stat = "mean")

  # At June: three steps back is March, May, June.
  expect_equal(unname(by_month(by_step, "SST_mean3")["6"]), mean(c(3, 5, 6)))
  # Three calendar months back is April, May, June, and April is absent.
  expect_equal(unname(by_month(by_cal, "SST_mean3month")["6"]), mean(c(5, 6)))
})

test_that("a calendar window can be short of observations and say so", {
  env <- rolling_field(c(1, 2, 3, 5, 6), months = c(1, 2, 3, 5, 6))

  strict <- rolling_covariate(env, "SST", n = 3, by = "month", stat = "mean",
                              min_obs = 3)

  # June's three-month window holds only two observations.
  expect_true(all(is.na(strict$SST_mean3month[strict$MONTH == 6])))
})

test_that("a year window spans twelve months across a year boundary", {
  env <- rolling_field(1:24, years = rep(c(2019L, 2020L), each = 12),
                       months = rep(1:12, 2))

  result <- rolling_covariate(env, "SST", n = 1, by = "year", stat = "mean")

  # December 2020 averages January to December 2020, values 13 to 24.
  last <- result[result$YEAR == 2020 & result$MONTH == 12, ]
  expect_equal(unique(last$SST_mean1year), mean(13:24))
})

test_that("each location is summarised over its own history", {
  grid <- expand.grid(x = c(-70, -69), y = 43)
  frames <- lapply(1:4, function(m) {
    f <- grid
    f$YEAR <- 2020L
    f$MONTH <- as.integer(m)
    f$DAY <- 1L
    f$SST <- ifelse(f$x < -69.5, m, 100 * m)
    f
  })
  env <- sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)

  result <- rolling_covariate(env, "SST", n = 2, stat = "mean")

  west <- sf::st_coordinates(env)[, 1] == -70
  expect_equal(result$SST_mean2[west & env$MONTH == 4], mean(c(3, 4)))
  expect_equal(result$SST_mean2[!west & env$MONTH == 4], mean(c(300, 400)))
})

test_that("missing values are skipped rather than poisoning the window", {
  env <- rolling_field(c(1, NA, 3, 4))

  result <- rolling_covariate(env, "SST", n = 3, stat = "mean")

  expect_equal(unname(by_month(result, "SST_mean3")[3]), mean(c(1, 3)))
})

test_that("a window with nothing in it is NA, not an infinite minimum", {
  # min() of an empty vector would be Inf with a warning, which would then look
  # like a real extreme downstream.
  env <- rolling_field(c(NA, NA, 3, 5))

  result <- rolling_covariate(env, "SST", n = 2, stat = c("min", "max"))

  expect_true(is.na(unname(by_month(result, "SST_min2")[1])))
  expect_true(is.na(unname(by_month(result, "SST_min2")[2])))
  expect_equal(unname(by_month(result, "SST_max2")[4]), 5)
})

test_that("several statistics and covariates come back in one call", {
  env <- rolling_field(c(1, 2, 3, 4))
  env$SSS <- 32 + env$MONTH

  result <- rolling_covariate(env, c("SST", "SSS"), n = 2,
                              stat = c("mean", "max"))

  expect_true(all(c("SST_mean2", "SST_max2", "SSS_mean2", "SSS_max2") %in%
                    names(result)))
})

test_that("suffixes can be given, one per statistic", {
  env <- rolling_field(c(1, 2, 3))

  result <- rolling_covariate(env, "SST", n = 2, stat = c("mean", "sd"),
                              suffix = c("_avg", "_var"))

  expect_true(all(c("SST_avg", "SST_var") %in% names(result)))
  expect_error(
    rolling_covariate(env, "SST", n = 2, stat = c("mean", "sd"),
                      suffix = "_only_one"),
    "one entry per statistic"
  )
})

test_that("an impossible window is rejected", {
  env <- rolling_field(c(1, 2, 3))
  expect_error(rolling_covariate(env, "SST", n = 0), "at least 1")
  expect_error(rolling_covariate(env, "SST", n = 2.5), "whole number")
})

test_that("a missing covariate is an error", {
  env <- rolling_field(c(1, 2, 3))
  expect_error(rolling_covariate(env, "NOPE"), "NOPE")
})
