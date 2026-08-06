test_that("cells on a front are at zero distance and others increase away from it", {
  # A step change in temperature at -69: a sharp north-south front down the
  # middle of the domain.
  env <- make_env(function(lon, lat, year, month) ifelse(lon > -69, 15, 10),
                  lon = seq(-70, -68, by = 0.1))

  result <- distance_to_front(env, "SST", quantile = 0.9)
  coords <- sf::st_coordinates(result)
  keep <- !is.na(result$SST_front_dist)
  found <- data.frame(lon = coords[keep, 1], distance = result$SST_front_dist[keep])

  expect_equal(min(found$distance), 0)

  # Distance must fall off monotonically with longitude as the front at -69 is
  # approached from the west.
  by_longitude <- aggregate(distance ~ lon, found, mean)
  by_longitude <- by_longitude[by_longitude$lon <= -69, ]
  by_longitude <- by_longitude[order(by_longitude$lon), ]
  expect_true(all(diff(by_longitude$distance) < 0))
})

test_that("distances are in kilometres by default and metres on request", {
  env <- make_env(function(lon, lat, year, month) ifelse(lon > -69, 15, 10),
                  lon = seq(-70, -68, by = 0.1))

  in_km <- distance_to_front(env, "SST")$SST_front_dist
  in_m <- distance_to_front(env, "SST", per = "m")$SST_front_dist

  # Same geometry, different unit. Compared as a ratio on the non-zero cells,
  # since 0 / 0 is not informative.
  nonzero <- which(!is.na(in_km) & in_km > 0)
  expect_equal(in_m[nonzero] / in_km[nonzero], rep(1000, length(nonzero)),
               tolerance = 1e-6)
})

test_that("an absolute threshold overrides the quantile", {
  # A ramp steepening toward the east, so the gradient varies continuously and
  # lowering the threshold genuinely admits more cells as frontal. A step field
  # would not test this: its gradient is zero everywhere except at the step, so
  # no threshold below that peak changes the answer.
  env <- make_env(function(lon, lat, year, month) (lon + 70)^2,
                  lon = seq(-70, -68, by = 0.1))

  strict <- distance_to_front(env, "SST", threshold = 0.03)
  permissive <- distance_to_front(env, "SST", threshold = 0.005)

  # More cells qualify as frontal, so nothing is further from one.
  expect_lt(mean(permissive$SST_front_dist, na.rm = TRUE),
            mean(strict$SST_front_dist, na.rm = TRUE))
})

test_that("a threshold above every gradient leaves no front to measure from", {
  env <- make_env(function(lon, lat, year, month) ifelse(lon > -69, 15, 10),
                  lon = seq(-70, -68, by = 0.1))

  result <- distance_to_front(env, "SST", threshold = 1e6)

  expect_true(all(is.na(result$SST_front_dist)))
})

test_that("record scope keeps one physical threshold across time steps", {
  # January has a sharp front; February is nearly flat. Under "record" the same
  # cutoff applies to both, so February - having no comparably sharp gradient -
  # correctly reports no front rather than inventing one.
  env <- make_env(
    function(lon, lat, year, month) {
      if (month == 1) ifelse(lon > -69, 15, 10) else 10 + (lon + 70) * 0.01
    },
    lon = seq(-70, -68, by = 0.1), months = 1:2
  )

  by_record <- distance_to_front(env, "SST", quantile = 0.95, scope = "record")
  february <- by_record$SST_front_dist[by_record$MONTH == 2]
  expect_true(all(is.na(february)))

  # Under "step" each month is thresholded against itself, so February gets
  # fronts by construction.
  by_step <- distance_to_front(env, "SST", quantile = 0.95, scope = "step")
  expect_false(all(is.na(by_step$SST_front_dist[by_step$MONTH == 2])))
})

test_that("a time step with no front is NA rather than zero", {
  # A completely uniform field has no gradient anywhere. Reporting 0 would say
  # "you are standing on a front", which is the opposite of the truth.
  env <- make_env(function(lon, lat, year, month) rep(10, length(lon)))

  expect_warning(result <- distance_to_front(env, "SST", threshold = 0.001),
                 "Spatially uniform")

  expect_true(all(is.na(result$SST_front_dist)))
})

test_that("the gradient helper column is not left behind", {
  env <- make_env(function(lon, lat, year, month) ifelse(lon > -69, 15, 10),
                  lon = seq(-70, -68, by = 0.1))

  result <- distance_to_front(env, "SST")

  expect_false(any(grepl("__front_grad", names(result))))
  expect_setequal(setdiff(names(result), names(env)), "SST_front_dist")
})

test_that("invalid thresholds are rejected", {
  env <- make_env()

  expect_error(distance_to_front(env, "SST", threshold = -1), "must be positive")
  expect_error(distance_to_front(env, "SST", quantile = 0), "between 0 and 1")
  expect_error(distance_to_front(env, "SST", quantile = 1), "between 0 and 1")
  expect_error(distance_to_front(env, "NOPE"), "Available")
})

test_that("distance to a contour is zero on it and grows away from it", {
  # Depth increasing linearly offshore, so the 100 m isobath is a straight line.
  env <- make_env(function(lon, lat, year, month) (lon + 70) * 100,
                  lon = seq(-70, -68, by = 0.1), var = "DEPTH")

  result <- distance_to_contour(env, "DEPTH", levels = 100)
  keep <- !is.na(result$DEPTH_dist_100)
  coords <- sf::st_coordinates(result)

  expect_true("DEPTH_dist_100" %in% names(result))
  expect_equal(min(result$DEPTH_dist_100, na.rm = TRUE), 0)

  # Depth = 100 at lon = -69, so distance must grow with separation from there.
  found <- data.frame(lon = coords[keep, 1], distance = result$DEPTH_dist_100[keep])
  by_longitude <- aggregate(distance ~ lon, found, mean)
  west <- by_longitude$distance[by_longitude$lon == min(by_longitude$lon)]
  middle <- by_longitude$distance[which.min(abs(by_longitude$lon + 69))]
  expect_gt(west, middle)
})

test_that("several isobaths produce one column each", {
  env <- make_env(function(lon, lat, year, month) (lon + 70) * 100,
                  lon = seq(-70, -68, by = 0.1), var = "DEPTH")

  result <- distance_to_isobath(env, levels = c(50, 100, 150))

  expect_true(all(c("isobath_dist_50", "isobath_dist_100", "isobath_dist_150")
                  %in% names(result)))
  # Each isobath sits somewhere different, so the distance fields differ.
  expect_false(isTRUE(all.equal(result$isobath_dist_50, result$isobath_dist_150)))
})

test_that("a level outside the covariate's range has no location", {
  env <- make_env(function(lon, lat, year, month) (lon + 70) * 100,
                  lon = seq(-70, -68, by = 0.1), var = "DEPTH")

  # Depth spans 0-200 m here, so there is no 5000 m isobath to measure to.
  result <- distance_to_isobath(env, levels = 5000)

  expect_true(all(is.na(result$isobath_dist_5000)))
})

test_that("the contour is found even when no cell sits exactly on it", {
  # Depth steps 0, 100, 200 with nothing at 150, so an exact-match test would
  # find no contour at all.
  env <- make_env(function(lon, lat, year, month) round((lon + 70) * 100, -2),
                  lon = seq(-70, -68, by = 0.1), var = "DEPTH")

  result <- distance_to_isobath(env, levels = 150)

  expect_false(all(is.na(result$isobath_dist_150)))
  expect_equal(min(result$isobath_dist_150, na.rm = TRUE), 0)
})

test_that("invalid levels are rejected", {
  env <- make_env(var = "DEPTH")

  expect_error(distance_to_contour(env, "DEPTH", levels = numeric(0)), "finite")
  expect_error(distance_to_contour(env, "DEPTH", levels = NA), "finite")
})
