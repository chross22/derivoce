# A Gulf of Maine grid, so distances can be checked against known geography
# rather than only against each other.
gom_env <- function(months = 1) {
  make_env(function(lon, lat, year, month) 10,
           lon = seq(-70, -66, by = 0.5), lat = seq(41, 44, by = 0.5),
           months = months)
}

test_that("distance to shore is positive and finite everywhere", {
  skip_if_not_installed("rnaturalearth")
  skip_if_not_installed("rnaturalearthdata")

  result <- distance_to_shore(gom_env())

  expect_true("shore_dist" %in% names(result))
  expect_false(any(is.na(result$shore_dist)))
  expect_true(all(result$shore_dist >= 0))
  expect_true(all(is.finite(result$shore_dist)))
})

test_that("offshore points are further from land than coastal ones", {
  skip_if_not_installed("rnaturalearth")
  skip_if_not_installed("rnaturalearthdata")

  result <- distance_to_shore(gom_env())
  coords <- sf::st_coordinates(result)

  # The southeastern corner is out past Georges Bank; the northwestern one is
  # against the Maine coast. If the two came out similar, the coastline was not
  # being found.
  offshore <- result$shore_dist[which.min((coords[, 1] + 66)^2 + (coords[, 2] - 41)^2)]
  coastal <- result$shore_dist[which.min((coords[, 1] + 70)^2 + (coords[, 2] - 44)^2)]

  expect_gt(offshore, coastal)
  expect_gt(offshore, 100)   # hundreds of km offshore
  expect_lt(coastal, 60)     # close to the Maine shore
})

test_that("distances are in kilometres, not metres or degrees", {
  skip_if_not_installed("rnaturalearth")
  skip_if_not_installed("rnaturalearthdata")

  result <- distance_to_shore(gom_env())

  # A 4-degree-wide shelf domain spans hundreds of km at most. Metres would put
  # these in the hundreds of thousands; degrees, below ten.
  expect_lt(max(result$shore_dist), 1000)
  expect_gt(max(result$shore_dist), 50)
})

test_that("the same location gets the same distance in every time step", {
  skip_if_not_installed("rnaturalearth")
  skip_if_not_installed("rnaturalearthdata")

  result <- distance_to_shore(gom_env(months = 1:3))

  first <- result$shore_dist[result$MONTH == 1]
  # Shore does not move, so a time-varying answer would mean the per-location
  # caching had mismatched its keys.
  expect_equal(result$shore_dist[result$MONTH == 2], first)
  expect_equal(result$shore_dist[result$MONTH == 3], first)
})

test_that("the coastline is padded beyond the data's extent", {
  skip_if_not_installed("rnaturalearth")
  skip_if_not_installed("rnaturalearthdata")

  # A box entirely offshore, whose nearest land lies outside it. Cropping the
  # coastline tightly to the data would find no land at all.
  offshore_only <- make_env(function(lon, lat, year, month) 10,
                            lon = seq(-67.5, -67, by = 0.25),
                            lat = seq(41, 41.5, by = 0.25))

  result <- distance_to_shore(offshore_only)

  expect_false(any(is.na(result$shore_dist)))
  expect_true(all(result$shore_dist > 0))
})

test_that("a margin too small to reach land is reported", {
  skip_if_not_installed("rnaturalearth")
  skip_if_not_installed("rnaturalearthdata")

  offshore_only <- make_env(function(lon, lat, year, month) 10,
                            lon = seq(-67.5, -67, by = 0.25),
                            lat = seq(41, 41.5, by = 0.25))

  expect_error(distance_to_shore(offshore_only, margin = 0.01), "No coastline")
})

test_that("the large coastline requires its data package", {
  skip_if_not_installed("rnaturalearth")
  skip_if(requireNamespace("rnaturalearthhires", quietly = TRUE),
          "rnaturalearthhires is installed, so the guard does not fire")

  expect_error(distance_to_shore(gom_env(), resolution = "large"),
               "rnaturalearthhires")
})
