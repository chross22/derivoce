# Decomposing a series into trend, seasonal cycle, and residual.
#
# The trend is fitted against *elapsed days*, not against step number, which is
# what lets a gap in the record leave a gap in the trend rather than compressing
# it. That has a consequence worth knowing when writing tests: calendar months
# are of unequal length, so a series that is a clean function of month number is
# not a clean function of time. A twelve-month sinusoid sampled on the first of
# each month is therefore not exactly orthogonal to a straight line in days, and
# a decomposition of it leaves a small real residual.
#
# So the fixtures below build series as functions of elapsed days where an exact
# answer is wanted, and the seasonal cases assert dominance rather than machine
# zero, with the reason stated at each one.

decompose_field <- function(value_fun, years = 2000:2009, months = 1:12,
                            lon = c(-70, -69.5), lat = 43) {
  grid <- expand.grid(x = lon, y = lat)
  stamps <- expand.grid(MONTH = months, YEAR = years)
  dates <- as.Date(paste(stamps$YEAR, stamps$MONTH, 1, sep = "-"))
  days <- as.numeric(dates - min(dates))

  frames <- lapply(seq_len(nrow(stamps)), function(i) {
    f <- grid
    f$YEAR <- as.integer(stamps$YEAR[i])
    f$MONTH <- as.integer(stamps$MONTH[i])
    f$DAY <- 1L
    f$SST <- value_fun(days[i], stamps$MONTH[i], f$x)
    f
  })
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}

cell_mean <- function(env, var) {
  key <- paste(round(sf::st_coordinates(env)[, 1], 6),
               round(sf::st_coordinates(env)[, 2], 6))
  ave <- tapply(env[[var]], key, mean, na.rm = TRUE)
  as.numeric(ave[match(key, names(ave))])
}


test_that("the parts add back up to the original", {
  env <- decompose_field(function(d, m, lon) {
    10 + 0.002 * d + 3 * sin(2 * pi * m / 12) + 0.4 * ((d * 7 + m * 13) %% 5)
  })

  out <- decompose_covariate(env, "SST")
  rebuilt <- cell_mean(out, "SST") + out$SST_trend + out$SST_seasonal +
    out$SST_residual

  expect_equal(rebuilt, out$SST)
})

test_that("a series linear in time is all trend, exactly", {
  env <- decompose_field(function(d, m, lon) 10 + 0.002 * d)

  out <- decompose_covariate(env, "SST")

  expect_true(all(abs(out$SST_residual) < 1e-8))
  expect_true(all(abs(out$SST_seasonal) < 1e-8))
  expect_gt(stats::sd(out$SST_trend), 1)
})

test_that("a repeating seasonal cycle is overwhelmingly seasonal", {
  # Not exactly, and the reason is the file header: a monthly sinusoid sampled
  # on unequal months is not quite orthogonal to a straight line in days.
  env <- decompose_field(function(d, m, lon) 10 + 5 * sin(2 * pi * m / 12))

  out <- decompose_covariate(env, "SST")

  expect_gt(stats::sd(out$SST_seasonal), 3)
  expect_lt(stats::sd(out$SST_trend), 0.05)
  expect_lt(stats::sd(out$SST_residual), 0.05)
  expect_equal(max(out$SST_seasonal), 5, tolerance = 0.02)
})

test_that("trend and seasonal are separated rather than confused", {
  env <- decompose_field(function(d, m, lon) {
    10 + 0.002 * d + 4 * sin(2 * pi * m / 12)
  })

  out <- decompose_covariate(env, "SST")

  # Both recovered, and almost nothing left over.
  expect_gt(stats::sd(out$SST_seasonal), 2)
  expect_gt(stats::sd(out$SST_trend), 0.5)
  expect_lt(stats::sd(out$SST_residual), 0.05)

  # The seasonal term repeats: every January carries the same value.
  january <- out$SST_seasonal[out$MONTH == 1]
  expect_equal(length(unique(round(january, 8))), 1L)
})

test_that("the trend is taken out before the seasonal cycle", {
  # A record starting mid-year and warming throughout gives the later months
  # more warm years than the earlier ones. Taking the seasonal cycle first
  # would bake that asymmetry into it as a step that is really the trend.
  env <- decompose_field(function(d, m, lon) 10 + 0.003 * d,
                         years = 2000:2004, months = 1:12)
  env <- env[!(env$YEAR == 2000 & env$MONTH < 7), ]

  out <- decompose_covariate(env, "SST")

  # A pure trend, so the seasonal term should be small despite the ragged start.
  expect_lt(stats::sd(out$SST_seasonal), 0.05 * stats::sd(out$SST_trend))
})

test_that("both centred parts have mean zero", {
  env <- decompose_field(function(d, m, lon) {
    10 + 0.002 * d + 2 * sin(2 * pi * m / 12)
  })

  out <- decompose_covariate(env, "SST")

  expect_equal(mean(out$SST_trend), 0, tolerance = 1e-8)
  expect_equal(mean(out$SST_seasonal), 0, tolerance = 1e-8)
})

test_that("slope is per year, not per day", {
  # Built to rise by exactly half a degree per 365.25 days.
  env <- decompose_field(function(d, m, lon) 10 + (0.5 / 365.25) * d)

  out <- decompose_covariate(env, "SST", components = "slope")

  expect_equal(unique(out$SST_slope), 0.5, tolerance = 1e-8)
})

test_that("slope sign follows the direction of change", {
  warming <- decompose_field(function(d, m, lon) 10 + 0.002 * d)
  cooling <- decompose_field(function(d, m, lon) 10 - 0.002 * d)

  expect_gt(unique(decompose_covariate(warming, "SST",
                                       components = "slope")$SST_slope), 0)
  expect_lt(unique(decompose_covariate(cooling, "SST",
                                       components = "slope")$SST_slope), 0)
})

test_that("each cell gets its own trend", {
  env <- decompose_field(function(d, m, lon) {
    10 + ifelse(lon < -69.7, 1, -1) * 0.002 * d
  })

  out <- decompose_covariate(env, "SST", components = "slope")

  west <- sf::st_coordinates(env)[, 1] < -69.7
  expect_gt(unique(out$SST_slope[west]), 0)
  expect_lt(unique(out$SST_slope[!west]), 0)
})

test_that("slope is refused for a curved trend", {
  env <- decompose_field(function(d, m, lon) 10 + 1e-6 * d^2)
  expect_error(
    decompose_covariate(env, "SST", degree = 2, components = "slope"),
    "degree = 1"
  )
})

test_that("a quadratic trend is captured by degree 2 and not by degree 1", {
  env <- decompose_field(function(d, m, lon) 10 + 1e-6 * d^2)

  linear <- decompose_covariate(env, "SST", components = c("trend", "residual"))
  curved <- decompose_covariate(env, "SST", degree = 2,
                                components = c("trend", "residual"))

  expect_gt(stats::sd(linear$SST_residual), stats::sd(curved$SST_residual))
  expect_true(all(abs(curved$SST_residual) < 1e-8))
})

test_that("a record too short for the trend is NA and says so", {
  env <- decompose_field(function(d, m, lon) 10 + 0.002 * d,
                         years = 2000L, months = 1L)

  expect_warning(out <- decompose_covariate(env, "SST", components = "trend"),
                 "fewest")
  expect_true(all(is.na(out$SST_trend)))
})

test_that("one year of monthly data cannot support a seasonal term", {
  env <- decompose_field(function(d, m, lon) 10 + m, years = 2000L)

  expect_warning(decompose_covariate(env, "SST"), "seasonal")
})

test_that("only the requested components are added", {
  env <- decompose_field(function(d, m, lon) 10 + 0.002 * d)

  out <- decompose_covariate(env, "SST", components = c("trend", "residual"))

  expect_true(all(c("SST_trend", "SST_residual") %in% names(out)))
  expect_false("SST_seasonal" %in% names(out))
  expect_false("SST_slope" %in% names(out))
})

test_that("missing values do not break the fit", {
  env <- decompose_field(function(d, m, lon) {
    10 + 0.002 * d + 2 * sin(2 * pi * m / 12)
  })
  env$SST[env$YEAR == 2005 & env$MONTH == 6] <- NA

  out <- decompose_covariate(env, "SST")

  expect_true(all(is.na(out$SST_residual[is.na(env$SST)])))
  expect_false(any(is.na(out$SST_trend)))
})

test_that("a gap in the record leaves a gap in the trend", {
  # Elapsed days rather than step positions: dropping three years must not
  # compress the trend into the years either side of the hole.
  full <- decompose_field(function(d, m, lon) 10 + 0.002 * d, months = 1L)
  gapped <- full[!(full$YEAR %in% 2003:2005), ]

  a <- decompose_covariate(full, "SST", components = "slope")
  b <- decompose_covariate(gapped, "SST", components = "slope")

  expect_equal(unique(a$SST_slope), unique(b$SST_slope), tolerance = 1e-8)
})

test_that("bad arguments are rejected", {
  env <- decompose_field(function(d, m, lon) 10 + 0.002 * d)
  expect_error(decompose_covariate(env, "SST", degree = 0), "at least 1")
  expect_error(decompose_covariate(env, "SST", degree = 1.5), "whole number")
  expect_error(decompose_covariate(env, "NOPE"), "NOPE")
})
