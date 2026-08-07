# Water-mass fractions, box anomalies, and the index dictionary.

flat_env <- function(temperature, salinity, months = 1) {
  env <- make_env(function(lon, lat, year, month) temperature, months = months)
  env$SSS <- salinity
  env
}

# ---- water mass fraction ---------------------------------------------------

slope_waters <- list(LSW = c(6, 34.4), WSW = c(12, 35.4))

test_that("a cell sitting on an endmember returns that endmember", {
  expect_true(all(water_mass_fraction(flat_env(6, 34.4), slope_waters)$lsw_frac == 1))
  expect_true(all(water_mass_fraction(flat_env(12, 35.4), slope_waters)$lsw_frac == 0))
})

test_that("a cell halfway along the mixing line returns one half", {
  # Halfway in both temperature and salinity, so the answer does not depend on
  # how the two axes are weighted against each other.
  result <- water_mass_fraction(flat_env(9, 34.9), slope_waters)

  expect_true(all(abs(result$lsw_frac - 0.5) < 1e-9))
})

test_that("the residual separates water on the mixing line from water off it", {
  on_line <- water_mass_fraction(flat_env(9, 34.9), slope_waters, residual = TRUE)
  off_line <- water_mass_fraction(flat_env(9, 33.0), slope_waters, residual = TRUE)

  expect_true(all(on_line$lsw_frac_residual < 1e-9))
  expect_gt(unique(off_line$lsw_frac_residual), 0.5)
})

test_that("fractions are clamped rather than extrapolated", {
  # Colder and fresher than either endmember. Reporting 1.4 would present an
  # endmember problem as a measurement.
  result <- water_mass_fraction(flat_env(0, 33), slope_waters)

  expect_true(all(result$lsw_frac >= 0 & result$lsw_frac <= 1))
})

test_that("degenerate endmembers are refused with a specific reason", {
  env <- flat_env(9, 34.9)

  expect_error(water_mass_fraction(env, list(A = c(6, 34), B = c(6, 34))),
               "identical")
  expect_error(water_mass_fraction(env, list(A = c(6, 34), B = c(6, 35))),
               "same temperature")
  expect_error(water_mass_fraction(env, list(A = c(6, 34), B = c(9, 34))),
               "same salinity")
  expect_error(water_mass_fraction(env, list(A = c(6, 34))), "exactly two")
})

test_that("the fraction column is named after the first endmember", {
  result <- water_mass_fraction(flat_env(9, 34.9), slope_waters)
  expect_true("lsw_frac" %in% names(result))

  named <- water_mass_fraction(flat_env(9, 34.9), slope_waters, name = "shelf")
  expect_true("shelf" %in% names(named))
})

# ---- box anomaly -----------------------------------------------------------

whole_box <- list(xmin = -70, xmax = -66, ymin = 41, ymax = 44)

test_that("with no reference the result is the box mean itself", {
  env <- make_env(function(lon, lat, year, month) month, months = 1:4)

  result <- box_anomaly(env, "SST", whole_box, reference = "none")

  expect_equal(result$SST_box, env$SST)
})

test_that("a record anomaly has zero mean over the time steps", {
  env <- make_env(function(lon, lat, year, month) month, months = 1:4)

  result <- box_anomaly(env, "SST", whole_box, reference = "record")

  expect_lt(abs(sum(unique(result$SST_box_anom))), 1e-9)
})

test_that("a climatology anomaly removes a repeating seasonal cycle exactly", {
  env <- make_env(function(lon, lat, year, month) month,
                  years = 2020:2022, months = 1:4)

  result <- box_anomaly(env, "SST", whole_box, reference = "climatology")

  expect_true(all(abs(result$SST_box_anom) < 1e-9))
})

test_that("the box restricts which cells contribute", {
  # Value depends on longitude, so a half-domain box must give a different mean
  # from the whole domain.
  env <- make_env(function(lon, lat, year, month) lon)

  whole <- box_anomaly(env, "SST", whole_box, reference = "none")$SST_box
  half <- box_anomaly(env, "SST", list(xmin = -70, xmax = -68, ymin = 41, ymax = 44),
                      reference = "none")$SST_box

  expect_false(isTRUE(all.equal(unique(whole), unique(half))))
})

test_that("an empty box errors and reports both extents", {
  env <- make_env()

  message <- tryCatch(box_anomaly(env, "SST", list(xmin = 0, xmax = 1,
                                                   ymin = 0, ymax = 1)),
                      error = conditionMessage)
  expect_match(message, "no grid points")
  expect_match(message, "data:")
})

test_that("a malformed box is refused", {
  env <- make_env()

  expect_error(box_anomaly(env, "SST", list(xmin = 0, xmax = 1)), "xmin")
  expect_error(box_anomaly(env, "SST", list(xmin = 1, xmax = 0, ymin = 0, ymax = 1)),
               "empty")
})

test_that("eastern_gom_salinity is box_anomaly over a fixed, inspectable box", {
  box <- eastern_gom_box()
  expect_named(box, c("xmin", "xmax", "ymin", "ymax"))

  env <- make_env(function(lon, lat, year, month) 32 + month,
                  lon = seq(-68, -66, by = 0.25), lat = seq(43, 44.5, by = 0.25),
                  months = 1:3)
  env$SSS <- env$SST

  result <- eastern_gom_salinity(env)
  expect_true("egom_salinity" %in% names(result))
  expect_equal(result$egom_salinity,
               box_anomaly(env, "SSS", box)$SSS_box_anom)
})

# ---- the dictionary --------------------------------------------------------

test_that("the dictionary lists every named index, with a source", {
  dictionary <- as.data.frame(derived_indices())

  expect_setequal(dictionary$name,
                  c("scotian_shelf_inflow", "northeast_channel_inflow",
                    "water_mass_fraction", "eastern_gom_salinity"))
  # Every entry must cite something and say which way it points.
  expect_true(all(nzchar(dictionary$source)))
  expect_true(all(nzchar(dictionary$sign)))
  expect_true(all(nzchar(dictionary$description)))
})

test_that("every documented index name is a function that exists", {
  # Catches a dictionary that drifts away from the code it describes.
  for (name in as.data.frame(derived_indices())$name) {
    expect_true(is.function(get(name, envir = asNamespace("derivoce"))),
                info = name)
  }
})

test_that("the dictionary prints without error", {
  expect_output(print(derived_indices()), "Regional indices")
  expect_output(print(derived_indices()), "scotian_shelf_inflow")
})

test_that("markdown rendering is a valid pipe table", {
  markdown <- derived_indices(markdown = TRUE)

  expect_type(markdown, "character")
  expect_length(markdown, 1)

  lines <- strsplit(markdown, "\n")[[1]]
  # Header, separator, and one row per index.
  expect_equal(length(lines), nrow(as.data.frame(derived_indices())) + 2)
  expect_true(all(startsWith(lines, "|")))
  expect_true(all(endsWith(lines, "|")))
  expect_match(lines[2], "---")
  # Every row must have the same number of cells as the header.
  cells <- vapply(lines, function(l) lengths(regmatches(l, gregexpr("|", l, fixed = TRUE))),
                  integer(1))
  expect_equal(length(unique(cells)), 1)
  expect_match(markdown, "`scotian_shelf_inflow()`", fixed = TRUE)
})
