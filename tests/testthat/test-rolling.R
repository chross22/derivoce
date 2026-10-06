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

  result <- rolling_covariate(env, "SST", n = 3, stat = "mean", by = "month")

  # March averages January, February and March.
  expect_equal(unname(by_month(result, "SST_mean3month")[3]), 2)
  expect_equal(unname(by_month(result, "SST_mean3month")[5]), 4)
})

test_that("early steps summarise the short window they have", {
  env <- rolling_field(c(1, 2, 3, 4, 5))

  result <- rolling_covariate(env, "SST", n = 3, stat = "mean", by = "month")

  expect_equal(unname(by_month(result, "SST_mean3month")[1]), 1)
  expect_equal(unname(by_month(result, "SST_mean3month")[2]), 1.5)
})

test_that("min_obs withholds a summary of too short a window", {
  env <- rolling_field(c(1, 2, 3, 4, 5))

  result <- rolling_covariate(env, "SST", n = 3, stat = "mean", min_obs = 3, by = "month")

  expect_true(all(is.na(result$SST_mean3[result$MONTH %in% 1:2])))
  expect_equal(unname(by_month(result, "SST_mean3month")[3]), 2)
})

test_that("each statistic computes what it says", {
  env <- rolling_field(c(1, 2, 3, 10, 5))

  result <- rolling_covariate(env, "SST", n = 3,
                              stat = c("mean", "min", "max", "sum", "median",
                                       "range"), by = "month")

  at4 <- function(col) unname(by_month(result, col)[4])
  expect_equal(at4("SST_mean3month"), mean(c(2, 3, 10)))
  expect_equal(at4("SST_min3month"), 2)
  expect_equal(at4("SST_max3month"), 10)
  expect_equal(at4("SST_sum3month"), 15)
  expect_equal(at4("SST_median3month"), 3)
  expect_equal(at4("SST_range3month"), 8)
})

test_that("a standard deviation needs two values", {
  env <- rolling_field(c(1, 2, 3, 4, 5))

  result <- rolling_covariate(env, "SST", n = 3, stat = "sd", by = "month")

  expect_true(is.na(unname(by_month(result, "SST_sd3month")[1])))
  expect_equal(unname(by_month(result, "SST_sd3month")[3]), stats::sd(1:3))
})

test_that("a window of one returns the value itself", {
  env <- rolling_field(c(4, 7, 1))

  result <- rolling_covariate(env, "SST", n = 1, stat = "mean", by = "month")

  expect_equal(result$SST_mean1month, result$SST)
})

test_that("a calendar window on a complete record is the last n months", {
  env <- rolling_field(c(1, 2, 3, 4, 5))

  result <- rolling_covariate(env, "SST", n = 3, by = "month", stat = "mean")

  expect_equal(unname(by_month(result, "SST_mean3month")["5"]), mean(c(3, 4, 5)))
})

test_that("a gap shortens a calendar window rather than reaching further back", {
  # Monthly series missing April. Three calendar months back from June is April,
  # May and June, and April is absent, so only two values are in the window. A
  # position count would have reached back to March and called it three months.
  env <- rolling_field(c(1, 2, 3, 5, 6), months = c(1, 2, 3, 5, 6))

  result <- rolling_covariate(env, "SST", n = 3, by = "month", stat = "mean")

  expect_equal(unname(by_month(result, "SST_mean3month")["6"]), mean(c(5, 6)))
})

test_that("counting positions in the record is refused, and a unit is required", {
  env <- rolling_field(c(1, 2, 3))

  expect_error(rolling_covariate(env, "SST", n = 2, by = "step"),
               "no longer counts positions")
  expect_error(rolling_covariate(env, "SST", n = 2), "needs `by`")
  expect_error(rolling_covariate(env, "SST", n = 2, by = "fortnight"),
               "must be one of")
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

  result <- rolling_covariate(env, "SST", n = 2, stat = "mean", by = "month")

  west <- sf::st_coordinates(env)[, 1] == -70
  expect_equal(result$SST_mean2[west & env$MONTH == 4], mean(c(3, 4)))
  expect_equal(result$SST_mean2[!west & env$MONTH == 4], mean(c(300, 400)))
})

test_that("missing values are skipped rather than poisoning the window", {
  env <- rolling_field(c(1, NA, 3, 4))

  result <- rolling_covariate(env, "SST", n = 3, stat = "mean", by = "month")

  expect_equal(unname(by_month(result, "SST_mean3month")[3]), mean(c(1, 3)))
})

test_that("a window with nothing in it is NA, not an infinite minimum", {
  # min() of an empty vector would be Inf with a warning, which would then look
  # like a real extreme downstream.
  env <- rolling_field(c(NA, NA, 3, 5))

  result <- rolling_covariate(env, "SST", n = 2, stat = c("min", "max"), by = "month")

  expect_true(is.na(unname(by_month(result, "SST_min2month")[1])))
  expect_true(is.na(unname(by_month(result, "SST_min2month")[2])))
  expect_equal(unname(by_month(result, "SST_max2month")[4]), 5)
})

test_that("several statistics and covariates come back in one call", {
  env <- rolling_field(c(1, 2, 3, 4))
  env$SSS <- 32 + env$MONTH

  result <- rolling_covariate(env, c("SST", "SSS"), n = 2,
                              stat = c("mean", "max"), by = "month")

  expect_true(all(c("SST_mean2month", "SST_max2month", "SSS_mean2month", "SSS_max2month") %in%
                    names(result)))
})

test_that("suffixes can be given, one per statistic", {
  env <- rolling_field(c(1, 2, 3))

  result <- rolling_covariate(env, "SST", n = 2, stat = c("mean", "sd"),
                              suffix = c("_avg", "_var"), by = "month")

  expect_true(all(c("SST_avg", "SST_var") %in% names(result)))
  expect_error(
    rolling_covariate(env, "SST", n = 2, stat = c("mean", "sd"),
                      suffix = "_only_one", by = "month"),
    "one entry per statistic"
  )
})

test_that("an impossible window is rejected", {
  env <- rolling_field(c(1, 2, 3))
  expect_error(rolling_covariate(env, "SST", n = 0, by = "month"), "at least 1")
  expect_error(rolling_covariate(env, "SST", n = 2.5, by = "month"), "whole number")
})

test_that("a missing covariate is an error", {
  env <- rolling_field(c(1, 2, 3))
  expect_error(rolling_covariate(env, "NOPE", by = "month"), "NOPE")
})

test_that("a calendar window on daily data never includes later days", {
  # Daily record over two months. The window is trailing, so on 1 January the
  # monthly maximum can only be that day's own value: counting by month alone
  # would sweep in the other thirty days of January, which are the future.
  days <- seq(as.Date("2020-01-01"), as.Date("2020-02-29"), by = "day")
  frame <- data.frame(
    x = -70, y = 43,
    YEAR = as.integer(format(days, "%Y")),
    MONTH = as.integer(format(days, "%m")),
    DAY = as.integer(format(days, "%d")),
    SST = seq_along(days)
  )
  env <- sf::st_as_sf(frame, coords = c("x", "y"), crs = 4326)

  result <- rolling_covariate(env, "SST", n = 1, by = "month", stat = "max")

  expect_equal(result$SST_max1month, result$SST[seq_along(days)][
    pmax(seq_along(days), 1)])
  expect_equal(result$SST_max1month[1], 1)
  expect_equal(result$SST_max1month[10], 10)
})
