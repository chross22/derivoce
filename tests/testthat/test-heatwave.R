# Marine heatwaves and cold spells.
#
# The threshold is a percentile of each cell's own climatology for that calendar
# month, so these tests need several years: with one year per month the
# climatology is a single value and nothing can exceed it.
#
# Note that a percentile threshold always has exceedances -- about a tenth of
# steps at the 90th percentile, by construction. An event is a *run* of them, so
# what min_steps does is the substance of the definition, not a detail.
#
# Two things are easy to get wrong and invisible if untested: a run must be
# consecutive in *time*, not merely in row order, so a cell missing from a step
# must break its event; and the cold-spell case is the same code with a sign
# flipped, exactly the reuse that silently reports every cold month as a
# heatwave.

# Deterministic, with a seasonal cycle and enough spread between years that a
# 90th percentile is a real threshold rather than the sample maximum.
baseline <- function(y, m) {
  10 + 2 * sin(2 * pi * m / 12) + 0.1 * ((y * 7 + m * 13) %% 11)
}

wave_field <- function(value_fun = baseline, years = 2000:2011, months = 1:12,
                       lon = c(-70, -69.5), lat = c(43, 43.5)) {
  grid <- expand.grid(x = lon, y = lat)
  frames <- list()
  for (y in years) {
    for (m in months) {
      f <- grid
      f$YEAR <- as.integer(y)
      f$MONTH <- as.integer(m)
      f$DAY <- 1L
      f$SST <- value_fun(y, m)
      frames[[length(frames) + 1]] <- f
    }
  }
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}

# Baseline everywhere, except the named year/month pairs which get `value`.
spike <- function(year, months, value, drop_to = NULL) {
  function(y, m) {
    if (y == year && m %in% months) value else baseline(y, m)
  }
}


test_that("a constant series has nothing above its own threshold", {
  env <- wave_field(function(y, m) 10)

  expect_warning(result <- marine_heatwave(env, "SST"), "[Ss]tatic")

  expect_false(any(result$mhw_event))
})

test_that("about a tenth of steps sit above the 90th percentile", {
  # Not a heatwave detector on its own: this is what the threshold means, and
  # why min_steps exists.
  env <- wave_field()

  result <- marine_heatwave(env, "SST", percentile = 0.9, min_steps = 1)

  expect_lt(mean(result$mhw_event), 0.25)
  expect_gt(mean(result$mhw_event), 0.02)
})

test_that("an injected extreme is flagged and is the most intense step", {
  env <- wave_field(spike(2011, 7, 30))

  result <- marine_heatwave(env, "SST", min_steps = 1)

  hot <- env$YEAR == 2011 & env$MONTH == 7
  expect_true(all(result$mhw_event[hot]))
  expect_equal(which.max(result$mhw_intensity), which(hot)[1])
})

test_that("intensity is measured from the climatology, not the threshold", {
  env <- wave_field(spike(2011, 7, 30))

  result <- marine_heatwave(env, "SST", min_steps = 1)

  hot <- which(env$YEAR == 2011 & env$MONTH == 7)[1]
  july <- env$MONTH == 7 & sf::st_coordinates(env)[, 1] == -70
  climatology <- mean(env$SST[july])
  expect_equal(result$mhw_intensity[hot], 30 - climatology)
})

test_that("a cold spell is the mirror image, not another heatwave", {
  env <- wave_field(spike(2011, 7, -5))

  result <- marine_heatwave(env, "SST", direction = "cold", min_steps = 1)

  cold <- env$YEAR == 2011 & env$MONTH == 7
  expect_true(all(result$mcs_event[cold]))
  # A departure below the climatology, so negative.
  expect_lt(result$mcs_intensity[which(cold)[1]], 0)
})

test_that("a hot step is not a cold spell", {
  env <- wave_field(spike(2011, 7, 30))
  result <- marine_heatwave(env, "SST", direction = "cold", min_steps = 1)
  hot <- env$YEAR == 2011 & env$MONTH == 7
  expect_false(any(result$mcs_event[hot]))
})

test_that("duration counts consecutive steps and repeats across the event", {
  env <- wave_field(spike(2011, 3:5, 30))

  result <- marine_heatwave(env, "SST", min_steps = 3)

  warm <- env$YEAR == 2011 & env$MONTH %in% 3:5
  expect_true(all(result$mhw_event[warm]))
  expect_gte(min(result$mhw_duration[warm]), 3L)
  # One event, so one duration across all its rows within a cell.
  expect_equal(length(unique(result$mhw_duration[warm])), 1L)
})

test_that("min_steps discards events that are too short", {
  env <- wave_field(spike(2011, 7, 30))

  kept <- marine_heatwave(env, "SST", min_steps = 1)
  dropped <- marine_heatwave(env, "SST", min_steps = 3)

  hot <- env$YEAR == 2011 & env$MONTH == 7
  expect_true(all(kept$mhw_event[hot]))
  expect_false(any(dropped$mhw_event[hot]))
})

test_that("max_gap joins events across a short interruption", {
  # Warm, warm, ordinary, warm, warm.
  env <- wave_field(spike(2011, c(3, 4, 6, 7), 30))

  split <- marine_heatwave(env, "SST", min_steps = 5, max_gap = 0)
  joined <- marine_heatwave(env, "SST", min_steps = 5, max_gap = 1)

  warm <- env$YEAR == 2011 & env$MONTH %in% c(3, 4, 6, 7)
  # Two runs of two cannot meet a five-step minimum separately.
  expect_false(any(split$mhw_event[warm]))
  # Bridged, they are one event of five including the ordinary month between.
  expect_true(all(joined$mhw_event[warm]))
  expect_gte(min(joined$mhw_duration[warm]), 5L)
})

test_that("a cell missing from a step breaks that cell's event", {
  env <- wave_field(spike(2011, c(3, 4, 6), 30))
  xy <- sf::st_coordinates(env)
  # One cell has no value for month 5; the other cells still do.
  env <- env[!(env$YEAR == 2011 & env$MONTH == 5 &
                 xy[, 1] == -70 & xy[, 2] == 43), ]

  result <- marine_heatwave(env, "SST", min_steps = 3)

  # For the gapped cell, months 3, 4 and 6 are adjacent in row order but not in
  # time. For every other cell they are two runs of 2 and 1. Neither reaches 3.
  warm <- env$YEAR == 2011 & env$MONTH %in% c(3, 4, 6)
  expect_false(any(result$mhw_event[warm]))
})

test_that("a step missing from the whole record closes the gap it left", {
  # Documented behaviour rather than desirable behaviour: the step is not in
  # the sequence at all, so its neighbours are adjacent and the event runs
  # through. Pinned here so it cannot change unnoticed.
  env <- wave_field(spike(2011, c(3, 4, 6), 30))
  env <- env[!(env$YEAR == 2011 & env$MONTH == 5), ]

  result <- marine_heatwave(env, "SST", min_steps = 3)

  warm <- env$YEAR == 2011 & env$MONTH %in% c(3, 4, 6)
  expect_true(all(result$mhw_event[warm]))
  expect_equal(unique(result$mhw_duration[warm]), 3L)
})

test_that("cumulative intensity is the total over the event", {
  env <- wave_field(spike(2011, 3:5, 30))

  result <- marine_heatwave(env, "SST", min_steps = 3)

  # One cell, not one longitude: -70 appears at two latitudes.
  xy <- sf::st_coordinates(env)
  cell <- xy[, 1] == -70 & xy[, 2] == 43
  warm <- env$YEAR == 2011 & env$MONTH %in% 3:5 & cell
  expect_equal(unique(result$mhw_cumulative[warm]),
               sum(result$mhw_intensity[warm]))
})

test_that("a long mild event and a short fierce one differ as expected", {
  mild <- wave_field(spike(2011, 3:8, 15))
  fierce <- wave_field(spike(2011, 3, 40))

  m <- marine_heatwave(mild, "SST", min_steps = 1)
  f <- marine_heatwave(fierce, "SST", min_steps = 1)

  expect_gt(max(m$mhw_duration, na.rm = TRUE), max(f$mhw_duration, na.rm = TRUE))
  expect_gt(max(f$mhw_intensity, na.rm = TRUE), max(m$mhw_intensity, na.rm = TRUE))
})

test_that("categories rise with how far past the threshold a value reaches", {
  env <- wave_field(spike(2011, 7, 40))

  result <- marine_heatwave(env, "SST", min_steps = 1)

  cats <- result$mhw_category[!is.na(result$mhw_category)]
  expect_true(all(cats >= 1L & cats <= 4L))
  hot <- which(env$YEAR == 2011 & env$MONTH == 7)[1]
  expect_equal(result$mhw_category[hot], max(cats))
})

test_that("each cell gets its own threshold", {
  # East runs 10 degrees warmer, and both have the same anomalous month.
  grid <- expand.grid(x = c(-70, -69), y = 43)
  frames <- list()
  for (y in 2000:2011) {
    f <- grid
    f$YEAR <- y
    f$MONTH <- 1L
    f$DAY <- 1L
    f$SST <- ifelse(f$x < -69.5, 5, 15) + 0.1 * ((y * 7) %% 11) +
      ifelse(y == 2011, 20, 0)
    frames[[length(frames) + 1]] <- f
  }
  env <- sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)

  result <- marine_heatwave(env, "SST", min_steps = 1)

  # Both cells flag the same year, despite very different means.
  hot <- env$YEAR == 2011
  expect_true(all(result$mhw_event[hot]))
})

test_that("a short baseline warns that the threshold is just the sample maximum", {
  env <- wave_field(years = 2000:2001, months = 1L)
  expect_warning(marine_heatwave(env, "SST"), "percentile")
})

test_that("only the requested measures are added", {
  env <- wave_field(spike(2011, 7, 30))

  result <- marine_heatwave(env, "SST", measures = c("event", "intensity"))

  expect_true(all(c("mhw_event", "mhw_intensity") %in% names(result)))
  expect_false("mhw_category" %in% names(result))
})

test_that("the prefix can be chosen", {
  env <- wave_field(spike(2011, 7, 30))
  result <- marine_heatwave(env, "SST", prefix = "hot")
  expect_true("hot_event" %in% names(result))
})

test_that("an impossible percentile is rejected", {
  env <- wave_field()
  expect_error(marine_heatwave(env, "SST", percentile = 1), "between 0 and 1")
  expect_error(marine_heatwave(env, "SST", percentile = 0), "between 0 and 1")
})

test_that("a missing covariate is an error", {
  env <- wave_field()
  expect_error(marine_heatwave(env, "NOPE"), "NOPE")
})
