# Every kind of column datamatch can attach, run through the whole package.
#
# datamatch serves four quite different things, and only the first behaves like
# an ordinary covariate:
#
#   * Copernicus variables, which vary in space and time;
#   * seafloor terrain from attach_bathymetry(), which is static -- identical in
#     every time step, so a lag reproduces it and a temporal gradient is zero;
#   * climate indices from attach_climate_index(), which are spatially uniform --
#     one basin-wide value per step, so a spatial gradient is zero everywhere;
#   * a `<var>_source` factor from fill_satellite_gaps(), which is not numeric
#     and cannot be differentiated at all.
#
# The last three are the ones that produce plausible-looking nonsense if a
# function does not notice what it has been handed. These tests assert the
# noticing, across the derivations added since the checks in test-checks.R.

datamatch_env <- function(years = 2000:2004) {
  lon <- seq(-70, -69, by = 0.25)
  lat <- seq(42.5, 43.5, by = 0.25)
  frames <- list()
  for (y in years) for (m in 1:12) {
    g <- expand.grid(x = lon, y = lat)
    g$YEAR <- as.integer(y); g$MONTH <- as.integer(m); g$DAY <- 1L
    # Copernicus: varies in space and time.
    g$SST <- 12 + 3 * sin(2 * pi * m / 12) + 0.1 * (y - 2000) + 2 * (g$x + 69.5)
    g$SSS <- 32 + 0.4 * (g$x + 69.5) + 0.05 * m
    g$UO <- 0.2 + 0.05 * (g$x + 69.5) + 0.01 * m
    g$VO <- 0.1 - 0.05 * (g$y - 43) + 0.01 * m
    # attach_bathymetry(): static.
    g$DEPTH <- 100 + (g$x + 70) * 20
    g$SLOPE <- 2 + (g$y - 43)
    # attach_climate_index(): spatially uniform.
    g$NAO <- 0.5 * m + 0.1 * (y - 2000)
    g$AMOC <- 17 + 0.05 * m
    # fill_satellite_gaps(): not numeric.
    g$SST_source <- factor(ifelse(m %% 2, "model", "satellite"))
    frames[[length(frames) + 1]] <- g
  }
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}

quiet <- function(expr) suppressWarnings(expr)


test_that("a non-numeric source column is skipped, not fumbled", {
  env <- datamatch_env()

  # vars = NULL sweeps every covariate; the factor must not stop it.
  out <- quiet(cell_anomaly(env))

  expect_true("SST_anom" %in% names(out))
  expect_false("SST_source_anom" %in% names(out))
})

test_that("naming the factor explicitly is an error that explains itself", {
  env <- datamatch_env(2000:2001)

  expect_error(cell_anomaly(env, "SST_source"), "not numeric")
  expect_error(rolling_covariate(env, "SST_source"), "not numeric")
  expect_error(decompose_covariate(env, "SST_source"), "not numeric")
})

test_that("static terrain warns in the new temporal derivations", {
  env <- datamatch_env(2000:2002)

  expect_warning(rolling_covariate(env, "DEPTH", n = 3), "[Ss]tatic")
  expect_warning(cell_anomaly(env, "DEPTH"), "[Ss]tatic")
  expect_warning(decompose_covariate(env, "DEPTH", components = "trend"),
                 "[Ss]tatic")
})

test_that("a spatially uniform index warns in the new spatial derivations", {
  env <- datamatch_env(2000:2001)

  expect_warning(front_frequency(env, "NAO"), "[Ss]patially uniform")
})

test_that("a climate index is still a legitimate temporal covariate", {
  # Uniform in space is not degenerate in time: lagging or decomposing an index
  # is exactly what it is for.
  env <- datamatch_env()

  expect_no_warning(out <- decompose_covariate(env, "NAO"))
  expect_false(any(is.na(out$NAO_trend)))
  expect_no_warning(rolling_covariate(env, "AMOC", n = 3, stat = "mean"))
})

test_that("index_series recovers the climate indices datamatch broadcast", {
  # The indices arrive as one value per step repeated on every row, which is
  # precisely what index_series() exists to undo.
  env <- datamatch_env(2000:2001)

  series <- index_series(env)

  expect_true(all(c("NAO", "AMOC") %in% names(series)))
  # Terrain varies across the grid, so it is not an index and must not appear.
  expect_false("DEPTH" %in% names(series))
  expect_false("SST" %in% names(series))
  expect_equal(nrow(series), 24)
})

test_that("naming a spatially varying column for index_series is refused", {
  env <- datamatch_env(2000:2001)
  expect_error(index_series(env, "SST"), "vary within a time step")
})

test_that("the physics derivations run on the Copernicus variables", {
  env <- datamatch_env(2000:2001)

  expect_no_error(quiet(flow_deformation(env)))
  expect_no_error(quiet(detect_eddies(env, measures = "in_eddy")))
  expect_no_error(quiet(potential_density(env)))
  expect_no_error(quiet(front_frequency(env, "SST")))
  expect_no_error(quiet(marine_heatwave(env, "SST")))
})


test_that("every datamatch variable has an EML unit", {
  known <- datamatch_units()

  expect_true(all(nzchar(known)))
  expect_equal(anyDuplicated(names(known)), 0L)
  # The four sources are all represented.
  expect_true(all(c("SST", "UO") %in% names(known)))        # Copernicus physics
  expect_true(all(c("CHL", "NO3") %in% names(known)))       # biogeochemistry
  expect_true(all(c("DEPTH", "TPI") %in% names(known)))     # bathymetry
  expect_true(all(c("NAO", "AMOC") %in% names(known)))      # climate indices
})

test_that("source units resolve without the caller supplying any", {
  env <- datamatch_env(2000:2001)
  env <- quiet(horizontal_gradient(env, "SST"))

  out <- eml_attributes(env, c("SST", "DEPTH", "NAO", "AMOC", "SST_grad"))
  unit <- function(v) out$unit[out$attributeName == v]

  expect_equal(unit("SST"), "celsius")
  expect_equal(unit("DEPTH"), "meter")
  expect_equal(unit("NAO"), "dimensionless")
  expect_equal(unit("AMOC"), "sverdrup")
  # And the derived unit is composed from the source one.
  expect_equal(unit("SST_grad"), "celsiusPerKilometer")
})

test_that("a supplied unit overrides the default for that column", {
  env <- datamatch_env(2000:2001)

  out <- eml_attributes(env, "SST", units = c(SST = "kelvin"))

  expect_equal(out$unit, "kelvin")
})

test_that("every datamatch unit is standard or declared custom", {
  # The invariant that keeps a document valid: EML rejects an undeclared unit.
  skip_if_not(requireNamespace("emld", quietly = TRUE) &&
                requireNamespace("xml2", quietly = TRUE), "emld not available")
  f <- system.file("tests/eml-2.1.1/eml-unitDictionary.xml", package = "emld")
  skip_if(!nzchar(f) || !file.exists(f), "unit dictionary not available")

  ids <- xml2::xml_attr(
    xml2::xml_find_all(xml2::read_xml(f), "//*[local-name()='unit']"), "id")
  ids <- ids[!is.na(ids)]

  undeclared <- setdiff(unname(datamatch_units()),
                        c(ids, eml_unit_registry()$id))
  expect_equal(undeclared, character(0))
})

test_that("the unit table has not fallen behind datamatch's catalogue", {
  # Hardcoded rather than read at runtime, because datamatch is deliberately
  # not a dependency. This is what catches the table going stale.
  skip_if_not(requireNamespace("datamatch", quietly = TRUE),
              "datamatch not installed")

  served <- unique(c(
    names(datamatch::copernicus_variables()),
    names(datamatch::ccmp_variables()),
    names(datamatch::hycom_variables()),
    names(datamatch::fvcom_variables()),
    unlist(lapply(datamatch::erddap_datasets(),
                  function(d) names(d$variables))),
    names(datamatch::bathymetry_variables()),
    names(datamatch::climate_indices())
  ))
  expect_equal(setdiff(served, names(datamatch_units())), character(0))
})

test_that("covariate_columns has not drifted from datamatch's again", {
  # The HOUR gap this file pins was exactly such a drift: the duplicate here
  # fell behind datamatch's own time_columns() and sub-daily fetches broke
  # silently. Comparing behaviour on a column zoo keeps the two in lockstep.
  skip_if_not(requireNamespace("datamatch", quietly = TRUE),
              "datamatch not installed")

  env <- datamatch_env(2000:2001)
  env$HOUR <- 12L
  env$SST_depth <- 0.5
  env$.datamatch_source <- "copernicus"

  expect_equal(derivoce:::covariate_columns(env),
               datamatch::covariate_columns(env))
})


# An unstructured mesh, as datamatch's accessFVCOM() returns: one row per mesh
# node, and the nodes irregularly spaced by design, because a mesh puts its
# resolution where the coastline is rather than on a lattice.
mesh_env <- function(years = 2000:2003) {
  set.seed(42)
  n <- 90
  nodes <- data.frame(x = runif(n, -70, -69), y = runif(n, 42.5, 43.5))
  frames <- list()
  for (y in years) for (m in 1:12) {
    g <- nodes
    g$YEAR <- as.integer(y); g$MONTH <- as.integer(m); g$DAY <- 1L
    g$SST <- 11 + 3 * sin(2 * pi * m / 12) + 0.1 * (y - 2000) + (g$x + 69.5)
    g$SSS <- 32 + 0.3 * (g$x + 69.5)
    g$UO <- 0.2 + 0.05 * (g$x + 69.5)
    g$VO <- 0.1 - 0.05 * (g$y - 43)
    frames[[length(frames) + 1]] <- g
  }
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}


test_that("every temporal derivation works on an unstructured mesh", {
  # These match points by coordinate and never build a raster, so a mesh is no
  # different to them from a lattice.
  mesh <- mesh_env()

  expect_no_error(lag_covariate(mesh, "SST", n = 1, by = "month"))
  expect_no_error(rolling_covariate(mesh, "SST", n = 3, stat = "mean"))
  expect_no_error(integrate_covariate(mesh, "SST"))
  expect_no_error(cell_anomaly(mesh, "SST"))
  expect_no_error(decompose_covariate(mesh, "SST"))
  expect_no_error(suppressWarnings(marine_heatwave(mesh, "SST")))
  expect_no_error(potential_density(mesh))
  expect_no_error(eml_attributes(mesh))
})

test_that("the derivations that need only coordinates work on a mesh too", {
  mesh <- mesh_env(2000:2001)

  expect_no_error(box_anomaly(mesh, "SST", box = list(
    xmin = -70, xmax = -69, ymin = 42.5, ymax = 43.5)))
})

test_that("a spatial derivation refuses a mesh, and says what to do", {
  mesh <- mesh_env(2000:2001)

  expect_error(horizontal_gradient(mesh, "SST"), "regular lon/lat grid")
  expect_error(flow_deformation(mesh), "regular lon/lat grid")
  expect_error(detect_eddies(mesh), "regular lon/lat grid")
  expect_error(distance_to_front(mesh, "SST"), "regular lon/lat grid")
  expect_error(ftle(mesh, integration_days = 3), "regular lon/lat grid")
})

test_that("the refusal names the mesh case and the way forward", {
  # Passed as a promise, rasterize_step()'s failure used to surface from inside
  # terra as "error in evaluating the argument 'x'", burying the explanation.
  mesh <- mesh_env(2000:2001)

  message <- tryCatch(horizontal_gradient(mesh, "SST"),
                      error = function(e) conditionMessage(e))

  expect_match(message, "unstructured", fixed = TRUE)
  expect_match(message, "accessFVCOM", fixed = TRUE)
  expect_match(message, "upscale_grid", fixed = TRUE)
  expect_false(grepl("evaluating the argument", message, fixed = TRUE))
})


# A sub-daily fetch, as accessCCMP(frequency = "6hourly") and
# accessHYCOM(frequency = "3hourly") return: an HOUR column beside the usual
# three. datamatch's own time_columns() includes HOUR, and this package's
# duplicate had drifted behind it, which collapsed the four hours of a day
# into one time step -- silently. Gradients rasterized whichever hour wrote
# last, lags drew from an arbitrary hour, and HOUR itself was swept up as a
# covariate. Everything below pins the repaired behaviour.
hourly_env <- function(days = 1:5, hours = c(0L, 6L, 12L, 18L)) {
  lon <- seq(-70, -69, by = 0.25); lat <- seq(42.5, 43.5, by = 0.25)
  frames <- list()
  for (d in days) for (h in hours) {
    g <- expand.grid(x = lon, y = lat)
    g$YEAR <- 2020L; g$MONTH <- 1L; g$DAY <- as.integer(d); g$HOUR <- h
    g$WSPD <- 8 + 2 * sin(2 * pi * h / 24) + 0.5 * d + 0.2 * (g$x + 69.5)
    g$UWND <- 5 + 0.1 * h + 0.1 * (g$x + 69.5)
    g$VWND <- 2 - 0.05 * h - 0.1 * (g$y - 43)
    frames[[length(frames) + 1]] <- g
  }
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}

wind_at <- function(env, d, h) {
  rows <- env$DAY == d & env$HOUR == h
  env[rows, ][order(sf::st_coordinates(env[rows, ])[, 1],
                    sf::st_coordinates(env[rows, ])[, 2]), ]
}


test_that("hours are distinct time steps, not duplicates of their day", {
  env <- hourly_env()

  steps <- derivoce:::time_steps(env)

  expect_equal(nrow(steps), 20)
  expect_true("HOUR" %in% names(steps))
  # Ordered within the day as well as across days.
  expect_equal(steps$HOUR[steps$DAY == 1], c(0L, 6L, 12L, 18L))
})

test_that("HOUR is a time column, not a covariate", {
  # datamatch's covariate_columns() excludes it; the duplicate here must agree,
  # or vars = NULL would lag and differentiate the clock.
  env <- hourly_env()

  expect_false("HOUR" %in% derivoce:::covariate_columns(env))
  # The input column survives, of course; what must not exist is an anomaly
  # of the clock.
  expect_false("HOUR_anom" %in% names(suppressWarnings(cell_anomaly(env))))
})

test_that("a one-step lag on hourly data is the previous hour", {
  env <- hourly_env()

  out <- lag_covariate(env, "WSPD", n = 1, by = "step")

  expect_equal(wind_at(out, 2, 6)$WSPD_lag1, wind_at(env, 2, 0)$WSPD)
  # And across the day boundary: hour 0 draws from the last hour of yesterday.
  expect_equal(wind_at(out, 2, 0)$WSPD_lag1, wind_at(env, 1, 18)$WSPD)
})

test_that("a calendar-day lag lands on the same hour of the previous day", {
  env <- hourly_env()

  out <- lag_covariate(env, "WSPD", n = 1, by = "day")

  # Hour 6 discriminates: the sinusoid gives hours 0 and 12 the same value,
  # so those two could not tell a same-hour match from a first-hour match.
  expect_equal(wind_at(out, 3, 6)$WSPD_lag1day, wind_at(env, 2, 6)$WSPD)
  expect_false(isTRUE(all.equal(wind_at(out, 3, 6)$WSPD_lag1day,
                                wind_at(env, 2, 0)$WSPD)))
})

test_that("spatial derivations see each hour's own field", {
  env <- hourly_env(days = 1, hours = c(0L, 12L))

  out <- horizontal_gradient(env, "UWND")

  # The zonal wind gradient is the same 0.1 per degree at both hours, but the
  # fields differ by 1.2 m/s -- so if hours collapsed, one of them rasterized
  # over the other and this comparison of the raw fields would fail.
  expect_equal(wind_at(out, 1, 12)$UWND - wind_at(out, 1, 0)$UWND,
               rep(1.2, 25))
  expect_false(any(is.na(wind_at(out, 1, 12)$UWND_grad[7])))
})

test_that("a rolling day-window on hourly data trails from the instant", {
  env <- hourly_env()

  out <- rolling_covariate(env, "WSPD", n = 1, by = "day", stat = "max")

  # The trailing 1-day window at day 3 hour 0 covers (day2 h6, day2 h12,
  # day2 h18, day3 h0], whose per-location maximum is day 2 hour 6. A window
  # computed on whole dates would sweep in the rest of day 3, and day 3 hour 6
  # is strictly larger everywhere, so the two disagree at every location.
  at <- function(d, h) wind_at(env, d, h)$WSPD
  expected <- pmax(at(2, 6), at(2, 12), at(2, 18), at(3, 0))

  expect_equal(wind_at(out, 3, 0)$WSPD_max1day, expected)
  expect_true(all(wind_at(out, 3, 0)$WSPD_max1day < at(3, 6)))
})

test_that("index_series keeps the hour", {
  env <- hourly_env(days = 1:2)
  env$NAO <- rep(seq_len(8), each = 25)  # one value per (day, hour) step

  series <- index_series(env, "NAO")

  expect_true("HOUR" %in% names(series))
  expect_equal(nrow(series), 8)
})

test_that("daily data is untouched by hour support", {
  # No HOUR column, no change: the canonical trio still keys everything.
  env <- datamatch_env(2000:2001)

  steps <- derivoce:::time_steps(env)
  expect_false("HOUR" %in% names(steps))
  expect_equal(nrow(steps), 24)
})
