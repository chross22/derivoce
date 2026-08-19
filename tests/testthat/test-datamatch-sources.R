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

  served <- c(names(datamatch::copernicus_variables()),
              names(datamatch::bathymetry_variables()),
              names(datamatch::climate_indices()))
  expect_equal(setdiff(served, names(datamatch_units())), character(0))
})
