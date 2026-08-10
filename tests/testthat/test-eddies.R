# Eddies as objects.
#
# The fixture is a velocity field with two counter-rotating vortices in known
# places, so polarity, position and count are all known by construction. The
# sign convention is the thing most worth pinning: cyclonic is counter-clockwise
# and positive vorticity, and getting it backwards would leave every test of
# "an eddy was found" passing while labelling every core the wrong way.

# A Gaussian vortex centred at (cx, cy). `spin` is +1 for counter-clockwise.
vortex <- function(lon, lat, cx, cy, spin, strength = 0.6, width = 0.35) {
  dx <- (lon - cx) * cos(cy * pi / 180)
  dy <- lat - cy
  decay <- exp(-(dx^2 + dy^2) / (2 * width^2))
  list(u = -spin * strength * dy * decay,
       v = spin * strength * dx * decay)
}

eddy_field <- function(spins = c(1, -1), centres = list(c(-69.5, 42.5),
                                                        c(-67.5, 43.5)),
                       months = 1:2,
                       lon = seq(-70.5, -66.5, by = 0.1),
                       lat = seq(41.8, 44.2, by = 0.1)) {
  frames <- lapply(months, function(m) {
    grid <- expand.grid(x = lon, y = lat)
    grid$UO <- 0
    grid$VO <- 0
    for (k in seq_along(spins)) {
      w <- vortex(grid$x, grid$y, centres[[k]][1], centres[[k]][2], spins[k])
      grid$UO <- grid$UO + w$u
      grid$VO <- grid$VO + w$v
    }
    grid$YEAR <- 2020L
    grid$MONTH <- as.integer(m)
    grid$DAY <- 1L
    grid
  })
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}

near <- function(env, cx, cy, tol = 0.25) {
  xy <- sf::st_coordinates(env)
  abs(xy[, 1] - cx) < tol & abs(xy[, 2] - cy) < tol
}


test_that("both vortex cores are found", {
  env <- eddy_field()

  out <- detect_eddies(env)

  expect_true(any(out$in_eddy[near(env, -69.5, 42.5)] == 1))
  expect_true(any(out$in_eddy[near(env, -67.5, 43.5)] == 1))
})

test_that("in_eddy is 1 inside and 0 outside, never NA", {
  env <- eddy_field()

  out <- detect_eddies(env, measures = "in_eddy")

  expect_true(all(out$in_eddy %in% c(0, 1)))
  # The far corner is quiet water, well away from either vortex.
  expect_equal(unique(out$in_eddy[near(env, -70.4, 41.9, tol = 0.15)]), 0)
})

test_that("polarity separates counter-clockwise from clockwise", {
  env <- eddy_field(spins = c(1, -1))

  out <- detect_eddies(env)

  ccw <- out$polarity[near(env, -69.5, 42.5) & out$in_eddy == 1]
  cw <- out$polarity[near(env, -67.5, 43.5) & out$in_eddy == 1]

  expect_true(all(ccw == 1))
  expect_true(all(cw == -1))
})

test_that("reversing the flow reverses every polarity", {
  forward <- detect_eddies(eddy_field(spins = c(1, -1)))
  reversed <- detect_eddies(eddy_field(spins = c(-1, 1)))

  a <- forward$polarity[!is.na(forward$polarity)]
  b <- reversed$polarity[!is.na(reversed$polarity)]
  expect_equal(sort(unique(a)), c(-1, 1))
  expect_equal(sort(unique(b)), c(-1, 1))
  # Same cells, opposite labels.
  both <- !is.na(forward$polarity) & !is.na(reversed$polarity)
  expect_true(mean(forward$polarity[both] == -reversed$polarity[both]) > 0.9)
})

test_that("cells outside an eddy have no eddy to describe", {
  env <- eddy_field()

  out <- detect_eddies(env)

  outside <- out$in_eddy == 0
  expect_true(all(is.na(out$polarity[outside])))
  expect_true(all(is.na(out$radius[outside])))
})

test_that("radius is a positive size in kilometres", {
  env <- eddy_field()

  out <- detect_eddies(env, per = "km")
  radii <- out$radius[!is.na(out$radius)]

  expect_true(all(radii > 0))
  # A vortex of this width is tens of kilometres across, not thousands.
  expect_true(all(radii > 5 & radii < 200))
})

test_that("a wider vortex gets a larger radius", {
  narrow <- eddy_field(spins = 1, centres = list(c(-68.5, 43)))
  wide <- narrow
  # Rebuild the same single vortex at twice the width.
  grid <- sf::st_coordinates(wide)
  w <- vortex(grid[, 1], grid[, 2], -68.5, 43, spin = 1, width = 0.7)
  wide$UO <- w$u
  wide$VO <- w$v

  a <- detect_eddies(narrow, measures = "radius")
  b <- detect_eddies(wide, measures = "radius")

  expect_gt(max(b$radius, na.rm = TRUE), max(a$radius, na.rm = TRUE))
})

test_that("min_cells discards small patches", {
  env <- eddy_field()

  loose <- detect_eddies(env, min_cells = 1, measures = "in_eddy")
  strict <- detect_eddies(env, min_cells = 400, measures = "in_eddy")

  expect_gt(sum(loose$in_eddy), sum(strict$in_eddy))
})

test_that("a flow with no rotation yields no eddies", {
  # Pure strain: stretched east, squeezed north, no spin anywhere.
  grid <- expand.grid(x = seq(-70, -68, by = 0.1), y = seq(42, 44, by = 0.1))
  grid$UO <- (grid$x + 69)
  grid$VO <- -(grid$y - 43)
  grid$YEAR <- 2020L; grid$MONTH <- 1L; grid$DAY <- 1L
  env <- sf::st_as_sf(grid, coords = c("x", "y"), crs = 4326)

  out <- detect_eddies(env, measures = "in_eddy")

  expect_equal(sum(out$in_eddy), 0)
})

test_that("distance to an eddy is zero inside it and rises with separation", {
  env <- eddy_field(spins = 1, centres = list(c(-68.5, 43)))

  out <- distance_to_eddy(env)

  inside <- out$eddy_dist == 0
  expect_true(any(inside, na.rm = TRUE))

  xy <- sf::st_coordinates(env)
  from_centre <- sqrt(((xy[, 1] + 68.5) * 80)^2 + ((xy[, 2] - 43) * 111)^2)
  keep <- !is.na(out$eddy_dist)
  expect_gt(stats::cor(from_centre[keep], out$eddy_dist[keep]), 0.8)
})

test_that("polarity narrows which eddies count as targets", {
  env <- eddy_field(spins = c(1, -1))

  cyc <- distance_to_eddy(env, polarity = "cyclonic")
  anti <- distance_to_eddy(env, polarity = "anticyclonic")

  expect_true("cyclonic_eddy_dist" %in% names(cyc))
  expect_true("anticyclonic_eddy_dist" %in% names(anti))

  # Zero distance is reached at the counter-clockwise core for one and at the
  # clockwise core for the other.
  ccw <- near(env, -69.5, 42.5)
  cw <- near(env, -67.5, 43.5)
  expect_equal(min(cyc$cyclonic_eddy_dist[ccw], na.rm = TRUE), 0)
  expect_gt(min(cyc$cyclonic_eddy_dist[cw], na.rm = TRUE), 0)
  expect_equal(min(anti$anticyclonic_eddy_dist[cw], na.rm = TRUE), 0)
})

test_that("a step with no qualifying eddy returns NA, not a distance to nothing", {
  grid <- expand.grid(x = seq(-70, -68, by = 0.1), y = seq(42, 44, by = 0.1))
  grid$UO <- (grid$x + 69)
  grid$VO <- -(grid$y - 43)
  grid$YEAR <- 2020L; grid$MONTH <- 1L; grid$DAY <- 1L
  env <- sf::st_as_sf(grid, coords = c("x", "y"), crs = 4326)

  out <- distance_to_eddy(env)

  expect_true(all(is.na(out$eddy_dist)))
})

test_that("metres and kilometres differ by the factor they should", {
  env <- eddy_field(spins = 1, centres = list(c(-68.5, 43)))

  km <- distance_to_eddy(env, per = "km")
  m <- distance_to_eddy(env, per = "m")
  keep <- !is.na(km$eddy_dist)

  expect_equal(m$eddy_dist[keep], km$eddy_dist[keep] * 1000)
})

test_that("only the requested measures are added", {
  env <- eddy_field()

  out <- detect_eddies(env, measures = c("in_eddy", "polarity"))

  expect_true(all(c("in_eddy", "polarity") %in% names(out)))
  expect_false("radius" %in% names(out))
})

test_that("a suffix renames the columns", {
  env <- eddy_field()
  out <- detect_eddies(env, measures = "in_eddy", suffix = "_ow")
  expect_true("in_eddy_ow" %in% names(out))
})

test_that("bad arguments are rejected", {
  env <- eddy_field()
  expect_error(detect_eddies(env, threshold = 0), "positive")
  expect_error(detect_eddies(env, threshold = -0.2), "positive")
  expect_error(detect_eddies(env, min_cells = 0), "at least 1")
  expect_error(detect_eddies(env, min_cells = 2.5), "whole number")
  env$VO <- NULL
  expect_error(detect_eddies(env), "VO")
})
