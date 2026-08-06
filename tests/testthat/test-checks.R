# A grid carrying the kinds of column datamatch now attaches alongside the
# Copernicus variables: static seafloor terrain, a basin-wide climate index, and
# the factor fill_satellite_gaps() uses to record where each value came from.
enriched_env <- function() {
  env <- make_env(function(lon, lat, year, month) lon + lat + month,
                  months = 1:3)
  coords <- sf::st_coordinates(env)
  # Static: depends on position only, identical in every time step.
  env$DEPTH <- 100 + (coords[, 1] + 70) * 20
  # Spatially uniform: one value per month, the same everywhere.
  env$NAO <- env$MONTH * 0.5
  env$SST_source <- factor(ifelse(env$MONTH == 1, "model", "satellite"))
  env
}

test_that("lagging a static covariate warns", {
  env <- enriched_env()

  expect_warning(result <- lag_covariate(env, "DEPTH"), "Static covariate")

  # The warning is the point, but the computation still runs and still returns
  # what it promises: the lag of a static column is that column, except in the
  # first step, which has no predecessor.
  later <- env$MONTH > 1
  expect_equal(result$DEPTH_lag1[later], env$DEPTH[later])
  expect_true(all(is.na(result$DEPTH_lag1[!later])))
})

test_that("integrating and differencing a static covariate warn too", {
  env <- enriched_env()

  expect_warning(integrate_covariate(env, "DEPTH"), "Static covariate")
  expect_warning(temporal_gradient(env, "DEPTH"), "Static covariate")
})

test_that("a temporal operation on a real covariate does not warn", {
  env <- enriched_env()

  expect_no_warning(lag_covariate(env, "SST"))
  expect_no_warning(integrate_covariate(env, "SST"))
})

test_that("a horizontal gradient of a spatially uniform covariate warns", {
  env <- enriched_env()

  expect_warning(result <- horizontal_gradient(env, "NAO"), "Spatially uniform")

  expect_true(all(is.na(result$NAO_grad) |
                    abs(result$NAO_grad) < 1e-9))
})

test_that("a spatial operation on a real covariate does not warn", {
  env <- enriched_env()

  expect_no_warning(horizontal_gradient(env, "SST"))
  # DEPTH is static through time but varies in space, which is exactly what a
  # horizontal gradient wants: this is slope, and warning here would be wrong.
  expect_no_warning(horizontal_gradient(env, "DEPTH"))
  expect_no_warning(distance_to_isobath(env, levels = 120))
})

test_that("the two degeneracies are checked independently", {
  env <- enriched_env()

  # Static in time, fine in space.
  expect_warning(lag_covariate(env, "DEPTH"), "Static covariate")
  expect_no_warning(horizontal_gradient(env, "DEPTH"))

  # Uniform in space, fine in time.
  expect_warning(horizontal_gradient(env, "NAO"), "Spatially uniform")
  expect_no_warning(lag_covariate(env, "NAO"))
})

test_that("a wrapper warns once, not once per internal call", {
  env <- enriched_env()

  # temporal_gradient() calls lag_covariate() internally and distance_to_front()
  # calls horizontal_gradient(); neither should double up.
  expect_warning(temporal_gradient(env, "DEPTH"), "Static covariate")
  expect_warning(distance_to_front(env, "NAO"), "Spatially uniform")
})

test_that("a single time step does not trigger the static warning", {
  # With one step every column is trivially constant through time. Warning here
  # would fire on every well-formed single-step object.
  env <- make_env(months = 1)

  expect_no_warning(lag_covariate(env, "SST"))
})

test_that("an explicitly named non-numeric column is an error", {
  env <- enriched_env()

  expect_error(horizontal_gradient(env, "SST_source"), "not numeric")
  expect_error(lag_covariate(env, "SST_source"), "not numeric")
})

test_that("vars = NULL skips non-numeric columns instead of failing", {
  env <- enriched_env()

  # The caller did not name SST_source, so sweeping it up and erroring would
  # make the NULL default unusable on any gap-filled object.
  result <- suppressWarnings(lag_covariate(env))

  expect_true("SST_lag1" %in% names(result))
  expect_false("SST_source_lag1" %in% names(result))
})

test_that("vars = NULL still warns about the degenerate columns it swept up", {
  env <- enriched_env()

  # This is the case the warning is really aimed at: the caller asked for
  # "everything" without knowing everything now includes seafloor terrain.
  expect_warning(lag_covariate(env), "Static covariate")
})

test_that("a missing column is still an error, and names the alternatives", {
  env <- enriched_env()

  expect_error(horizontal_gradient(env, "NOPE"), "not present")
  expect_error(horizontal_gradient(env, "NOPE"), "SST")
})

test_that("an all-NA column is not reported as degenerate", {
  # No evidence either way, and a column of NAs has a more obvious problem than
  # the one this warning is about.
  env <- make_env(months = 1:3)
  env$EMPTY <- NA_real_

  expect_no_warning(lag_covariate(env, "EMPTY"))
  expect_no_warning(horizontal_gradient(env, "EMPTY"))
})
