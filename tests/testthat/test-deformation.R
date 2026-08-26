# Velocity gradient diagnostics.
#
# The signs are the easy thing to get wrong: vorticity and shear strain differ
# only in the sign between the same two derivatives, and terra's rows run north
# first, so a flipped y kernel would leave every identity below intact while
# reversing the rotation. The tests therefore check named flow patterns whose
# answer is known by inspection, not just the algebra.

deformation_field <- function(u_fun, v_fun, months = 1:3,
                              lon = seq(-70, -68, by = 0.25),
                              lat = seq(42, 44, by = 0.25)) {
  frames <- lapply(months, function(m) {
    grid <- expand.grid(x = lon, y = lat)
    grid$UO <- u_fun(grid$x, grid$y)
    grid$VO <- v_fun(grid$x, grid$y)
    grid$YEAR <- 2020
    grid$MONTH <- m
    grid$DAY <- 1L
    grid
  })
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}

zero <- function(lon, lat) rep(0, length(lon))
interior <- function(x) x[is.finite(x)]

all_measures <- c("vorticity", "divergence", "normal_strain", "shear_strain",
                  "strain_rate", "okubo_weiss")


test_that("a uniform current is not deforming, and says so", {
  env <- deformation_field(function(lon, lat) rep(0.4, length(lon)),
                           function(lon, lat) rep(-0.2, length(lon)))

  # A spatially uniform field has no gradient to take, which the package warns
  # about rather than silently returning a column of zeroes.
  expect_warning(result <- flow_deformation(env, measures = all_measures),
                 "uniform")

  for (m in all_measures) {
    expect_true(all(abs(interior(result[[m]])) < 1e-12),
                info = paste(m, "should vanish in a uniform flow"))
  }
})

test_that("counter-clockwise rotation has positive vorticity", {
  # u = -(y - y0), v = +(x - x0) turns anticlockwise seen from above.
  env <- deformation_field(function(lon, lat) -(lat - 43),
                           function(lon, lat) (lon + 69))

  result <- flow_deformation(env, measures = all_measures)

  expect_true(all(interior(result$vorticity) > 0))
  # Rotation without deformation: the strains cancel and Okubo-Weiss is negative.
  expect_true(all(abs(interior(result$strain_rate)) <
                    abs(interior(result$vorticity))))
  expect_true(all(interior(result$okubo_weiss) < 0))
})

test_that("clockwise rotation reverses the sign", {
  env <- deformation_field(function(lon, lat) (lat - 43),
                           function(lon, lat) -(lon + 69))

  result <- flow_deformation(env)

  expect_true(all(interior(result$vorticity) < 0))
  expect_true(all(interior(result$okubo_weiss) < 0))
})

test_that("outflow is positive divergence and inflow negative", {
  out <- deformation_field(function(lon, lat) (lon + 69),
                           function(lon, lat) (lat - 43))
  into <- deformation_field(function(lon, lat) -(lon + 69),
                            function(lon, lat) -(lat - 43))

  expect_true(all(interior(flow_deformation(out)$divergence) > 0))
  expect_true(all(interior(flow_deformation(into)$divergence) < 0))
})

test_that("a pure strain field deforms without rotating", {
  # Stretched east, squeezed north: no spin, so Okubo-Weiss is positive.
  env <- deformation_field(function(lon, lat) (lon + 69),
                           function(lon, lat) -(lat - 43))

  result <- flow_deformation(env, measures = all_measures)

  expect_true(all(abs(interior(result$vorticity)) < 1e-12))
  expect_true(all(interior(result$okubo_weiss) > 0))
  expect_true(all(interior(result$strain_rate) > 0))
})

test_that("Okubo-Weiss is the strain squared minus the vorticity squared", {
  env <- deformation_field(function(lon, lat) (lon + 69) * (lat - 43),
                           function(lon, lat) (lon + 69)^2)

  result <- flow_deformation(env, measures = all_measures)
  keep <- is.finite(result$okubo_weiss)

  expect_equal(
    result$okubo_weiss[keep],
    result$normal_strain[keep]^2 + result$shear_strain[keep]^2 -
      result$vorticity[keep]^2
  )
})

test_that("strain rate is the magnitude of the two strain components", {
  env <- deformation_field(function(lon, lat) (lon + 69) * (lat - 43),
                           function(lon, lat) (lat - 43)^2)

  result <- flow_deformation(env, measures = all_measures)
  keep <- is.finite(result$strain_rate)

  expect_equal(
    result$strain_rate[keep],
    sqrt(result$normal_strain[keep]^2 + result$shear_strain[keep]^2)
  )
})

test_that("derivatives are per second, so the magnitude is physical", {
  # v increases by 1 m/s per degree of longitude. At 43N a degree is
  # 111320 * cos(43) metres, so dv/dx is the reciprocal of that. u varies with
  # longitude only, so du/dy is zero and vorticity is dv/dx alone.
  env <- deformation_field(function(lon, lat) 0.1 * (lon + 69),
                           function(lon, lat) (lon + 69))

  result <- flow_deformation(env, measures = "vorticity")
  env$lat <- sf::st_coordinates(env)[, 2]
  at_43 <- result$vorticity[abs(env$lat - 43) < 1e-9]
  at_43 <- at_43[is.finite(at_43)]

  expected <- 1 / (111320 * cos(43 * pi / 180))
  expect_equal(at_43, rep(expected, length(at_43)), tolerance = 1e-6)
})

test_that("the Rossby number is vorticity scaled by the Coriolis parameter", {
  env <- deformation_field(function(lon, lat) -(lat - 43),
                           function(lon, lat) (lon + 69))

  result <- flow_deformation(env, measures = c("vorticity", "rossby"))
  keep <- is.finite(result$rossby)

  omega <- 7.2921e-5
  f <- 2 * omega * sin(sf::st_coordinates(env)[keep, 2] * pi / 180)
  expect_equal(result$rossby[keep], result$vorticity[keep] / f)

  # Northern hemisphere: f is positive, so the sign is carried through.
  expect_true(all(result$rossby[keep] > 0))
})

test_that("only the requested measures are added", {
  env <- deformation_field(function(lon, lat) (lon + 69),
                           function(lon, lat) (lat - 43))

  result <- flow_deformation(env, measures = c("divergence", "vorticity"))

  expect_true(all(c("divergence", "vorticity") %in% names(result)))
  expect_false("okubo_weiss" %in% names(result))
  expect_false("strain_rate" %in% names(result))
})

test_that("a suffix renames the columns", {
  env <- deformation_field(function(lon, lat) (lon + 69),
                           function(lon, lat) (lat - 43))

  result <- flow_deformation(env, measures = "divergence", suffix = "_srf")

  expect_true("divergence_srf" %in% names(result))
  expect_false("divergence" %in% names(result))
})

test_that("an unknown measure is rejected", {
  env <- deformation_field(zero, zero)
  expect_error(flow_deformation(env, measures = "enstrophy"))
})

test_that("a missing velocity column is an error naming what is available", {
  env <- deformation_field(zero, zero)
  env$VO <- NULL
  expect_error(flow_deformation(env), "VO")
})

test_that("the outermost ring is NA, as for every central difference", {
  env <- deformation_field(function(lon, lat) (lon + 69),
                           function(lon, lat) (lat - 43))

  result <- flow_deformation(env, measures = "divergence")

  expect_true(any(is.na(result$divergence)))
  expect_true(any(is.finite(result$divergence)))
})


test_that("a projected grid only stops the measure that needs latitude", {
  # The Rossby number is the one measure that cannot be computed without
  # latitude. Building all seven layers and then subsetting would evaluate it
  # on every call, so a request for vorticity alone on a projected grid would
  # fail citing a measure the caller never asked for.
  grid <- expand.grid(x = seq(0, 4e4, by = 1e4), y = seq(0, 4e4, by = 1e4))
  grid$UO <- grid$y / 1e5
  grid$VO <- -grid$x / 1e5
  grid$YEAR <- 2020
  grid$MONTH <- 1
  grid$DAY <- 1L
  env <- sf::st_as_sf(grid, coords = c("x", "y"), crs = "EPSG:32619")

  expect_no_error(flow_deformation(env, measures = all_measures))
  expect_error(flow_deformation(env, measures = "rossby"), "projected grid")
})
