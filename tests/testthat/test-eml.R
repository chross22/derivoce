# EML metadata for derived columns.
#
# The value of this is entirely in being right, so the tests are about honesty
# rather than coverage: a unit must be either correct or absent, never guessed;
# a column this package did not create must say so; and any unit outside EML's
# standard dictionary must come back in the custom declarations, because a
# document using an undeclared unit fails validation.

eml_field <- function() {
  # Large enough for a central difference: the spatial functions refuse a grid
  # with fewer than two distinct longitudes or latitudes, and rightly.
  grid <- expand.grid(x = seq(-70, -69, by = 0.25), y = seq(42.5, 43.5, by = 0.25))
  frames <- lapply(1:3, function(m) {
    f <- grid
    f$YEAR <- 2020L
    f$MONTH <- as.integer(m)
    f$DAY <- 1L
    f$SST <- 10 + m + f$x + 0.5 * f$y
    f$UO <- 0.2 + 0.05 * (f$x + 69.5)
    f$VO <- 0.1 - 0.05 * (f$y - 43)
    # Not a variable datamatch serves, so nothing is known about its units.
    f$TEMP_INSITU <- 9 + m + f$x
    f
  })
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}

# The 195 unit ids EML 2.1.1 recognises, if emld is around to check against.
standard_units <- function() {
  if (!requireNamespace("emld", quietly = TRUE) ||
      !requireNamespace("xml2", quietly = TRUE)) {
    return(NULL)
  }
  f <- system.file("tests/eml-2.1.1/eml-unitDictionary.xml", package = "emld")
  if (!nzchar(f) || !file.exists(f)) return(NULL)
  ids <- xml2::xml_attr(
    xml2::xml_find_all(xml2::read_xml(f), "//*[local-name()='unit']"), "id")
  ids[!is.na(ids)]
}


test_that("the table has the columns EML::set_attributes expects", {
  env <- eml_field()
  out <- eml_attributes(env, "SST")

  expect_true(all(c("attributeName", "attributeDefinition", "measurementScale",
                    "domain", "unit", "numberType") %in% names(out)))
  expect_equal(out$attributeName, "SST")
})

test_that("a derived column is described by what it is, not by its name", {
  env <- eml_field()
  env <- eke(env)

  out <- eml_attributes(env, "EKE")

  expect_match(out$attributeDefinition, "[Ee]ddy kinetic energy")
  expect_equal(out$unit, "metersSquaredPerSecondSquared")
})

test_that("a gradient's unit is built from the source column's unit", {
  env <- eml_field()
  env <- horizontal_gradient(env, "SST")

  out <- eml_attributes(env, "SST_grad", units = c(SST = "celsius"))

  expect_equal(out$unit, "celsiusPerKilometer")
  expect_match(out$attributeDefinition, "SST")
})

test_that("an unknown source unit gives no unit rather than a guess", {
  env <- eml_field()
  env <- horizontal_gradient(env, "TEMP_INSITU")

  out <- eml_attributes(env, "TEMP_INSITU_grad")

  expect_true(is.na(out$unit))
  # But it still says what the quantity is.
  expect_match(out$attributeDefinition, "gradient")
})

test_that("a column this package did not create says so", {
  env <- eml_field()

  out <- eml_attributes(env, "TEMP_INSITU")

  expect_match(out$attributeDefinition, "not derived by derivoce")
  expect_true(is.na(out$unit))
})

test_that("a datamatch column is still not derived here, but its unit is known", {
  # The two are independent: derivoce did not create SST, and yet it knows what
  # SST holds, because datamatch serves it.
  env <- eml_field()

  out <- eml_attributes(env, "SST")

  expect_match(out$attributeDefinition, "not derived by derivoce")
  expect_equal(out$unit, "celsius")
})

test_that("a supplied unit is used for the source column itself", {
  env <- eml_field()
  out <- eml_attributes(env, "SST", units = c(SST = "celsius"))
  expect_equal(out$unit, "celsius")
})

test_that("suffixes are matched longest-first, not greedily", {
  env <- eml_field()
  env <- distance_to_front(env, "SST")
  env <- rolling_covariate(env, "SST", n = 2, stat = "mean")

  front <- eml_attributes(env, "SST_front_dist")
  roll <- eml_attributes(env, "SST_mean2", units = c(SST = "celsius"))

  expect_equal(front$unit, "kilometer")
  expect_match(front$attributeDefinition, "front")
  expect_equal(roll$unit, "celsius")
  expect_match(roll$attributeDefinition, "trailing window")
})

test_that("a lag keeps the source unit and names its source", {
  env <- eml_field()
  env <- lag_covariate(env, "SST", n = 1, by = "month")

  out <- eml_attributes(env, "SST_lag1month", units = c(SST = "celsius"))

  expect_equal(out$unit, "celsius")
  expect_match(out$attributeDefinition, "lagged")
})

test_that("a z-score is dimensionless whatever its source was", {
  env <- eml_field()
  env <- cell_anomaly(env, "SST", reference = "record", standardize = TRUE)

  out <- eml_attributes(env, "SST_z", units = c(SST = "celsius"))

  expect_equal(out$unit, "dimensionless")
})

test_that("flow diagnostics carry the units they are computed in", {
  env <- eml_field()
  env <- flow_deformation(env, measures = c("vorticity", "okubo_weiss"))

  out <- eml_attributes(env, c("vorticity", "okubo_weiss"))

  expect_equal(out$unit[out$attributeName == "vorticity"], "perSecond")
  expect_equal(out$unit[out$attributeName == "okubo_weiss"], "perSecondSquared")
})

test_that("every column is described when none are named", {
  env <- eml_field()
  env <- current_speed(env)

  out <- eml_attributes(env)

  expect_true(all(c("SST", "UO", "VO", "speed") %in% out$attributeName))
  expect_equal(out$unit[out$attributeName == "speed"], "metersPerSecond")
})


test_that("custom units are returned for exactly the non-standard ones used", {
  env <- eml_field()
  env <- flow_deformation(env, measures = "vorticity")
  env <- current_speed(env)

  attributes <- eml_attributes(env)
  custom <- eml_custom_units(attributes)

  # perSecond is not in EML's dictionary and must be declared.
  expect_true("perSecond" %in% custom$id)
  # metersPerSecond is, so it must not be.
  expect_false("metersPerSecond" %in% custom$id)
})

test_that("no custom units are needed when every unit is standard", {
  env <- eml_field()
  env <- current_speed(env)

  custom <- eml_custom_units(eml_attributes(env, "speed"))

  expect_equal(nrow(custom), 0)
})

test_that("custom declarations carry what EML requires of them", {
  custom <- eml_custom_units("perSecond")

  expect_equal(nrow(custom), 1)
  expect_true(all(c("id", "unitType", "parentSI", "multiplierToSI",
                    "description") %in% names(custom)))
  expect_true(nzchar(custom$description))
})

test_that("every unit used is either standard or declared as custom", {
  # The invariant the whole thing rests on: a document using an undeclared
  # non-standard unit does not validate.
  ids <- standard_units()
  skip_if(is.null(ids), "emld unit dictionary not available")

  declared <- eml_unit_registry()$id
  registry_units <- unlist(lapply(eml_fixed_registry(), function(e) e$unit))
  registry_units <- unique(registry_units[!is.na(registry_units)])

  undeclared <- setdiff(registry_units, c(ids, declared))
  expect_equal(undeclared, character(0))
})

test_that("declared custom units do not duplicate standard ones", {
  ids <- standard_units()
  skip_if(is.null(ids), "emld unit dictionary not available")

  expect_equal(intersect(eml_unit_registry()$id, ids), character(0))
})

test_that("column classes line up with the attributes", {
  env <- eml_field()
  env <- current_speed(env)

  attributes <- eml_attributes(env)
  classes <- eml_col_classes(env)

  expect_equal(length(classes), nrow(attributes))
  expect_true(all(classes %in% c("numeric", "character", "factor", "Date")))
})

test_that("a logical column is numeric to EML, not character", {
  env <- eml_field()
  env$flag <- TRUE

  expect_equal(eml_col_classes(env, "flag"), "numeric")
  expect_equal(eml_attributes(env, "flag")$domain, "numericDomain")
})

test_that("a missing column is an error naming what is available", {
  env <- eml_field()
  expect_error(eml_attributes(env, "NOPE"), "NOPE")
})

# covariate_columns() is kept identical to datamatch's, since the two describe
# the same thing about the same objects and a reader should not have to check
# which package's version they are looking at.

test_that("covariate_columns matches datamatch: source in, depth out", {
  d <- sf::st_as_sf(
    data.frame(x = c(1, 2), y = c(1, 2), SST = c(4, 5),
               SST_source = c("satellite", "model"), BOTS_depth = c(80, 90),
               YEAR = 2010L, MONTH = 6L, DAY = 1L),
    coords = c("x", "y"), crs = 4326)

  cols <- covariate_columns(d)
  # A source tag travels with the variable it describes.
  expect_true("SST_source" %in% cols)
  # A model level is not a measurement, and its mean is not a depth.
  expect_false("BOTS_depth" %in% cols)
  expect_false(".datamatch_source" %in% cols)
})

test_that("DEPTH is a covariate, not provenance", {
  # distance_to_isobath() works off a DEPTH column from attach_bathymetry().
  # The bookkeeping pattern needs the underscore, so this must survive.
  d <- sf::st_as_sf(
    data.frame(x = c(1, 2), y = c(1, 2), DEPTH = c(100, 120),
               YEAR = 2010L, MONTH = 6L, DAY = 1L),
    coords = c("x", "y"), crs = 4326)
  expect_true("DEPTH" %in% covariate_columns(d))
})

test_that("an archive still describes a depth column", {
  # Narrowing covariate_columns() must not narrow what gets deposited:
  # datamatch::write_eml() documents <var>_depth with metres as its unit, and a
  # column in a deposited table with nothing describing it is worse than one
  # that is merely not a covariate.
  d <- sf::st_as_sf(
    data.frame(x = c(1, 2), y = c(1, 2), SST = c(4, 5), BOTS_depth = c(80, 90),
               YEAR = 2010L, MONTH = 6L, DAY = 1L),
    coords = c("x", "y"), crs = 4326)

  expect_true("BOTS_depth" %in% documented_columns(d))
  expect_true("BOTS_depth" %in% eml_attributes(d)$attributeName)
  # col_classes is positional, to line up with the attribute table EML wants,
  # so the check is that it covers the same columns rather than names them.
  expect_length(eml_col_classes(d), length(documented_columns(d)))
  expect_equal(eml_col_classes(d)[match("BOTS_depth", documented_columns(d))],
               "numeric")
})
