# Potential density.
#
# The coefficients are the whole function, and a mistyped digit would shift
# every value by a small constant that still looks like a plausible density.
# The published one-atmosphere check values are therefore the primary test,
# with the physics (a freshwater maximum near 4 C, saltier being denser) as the
# check that the polynomial has the right shape rather than merely the right
# value at four points.

ts_field <- function(t, s, months = seq_along(t)) {
  grid <- expand.grid(x = c(-70, -69.5), y = c(43, 43.5))
  frames <- lapply(seq_along(t), function(i) {
    f <- grid
    f$YEAR <- 2020L
    f$MONTH <- as.integer(months[i])
    f$DAY <- 1L
    f$SST <- t[i]
    f$SSS <- s[i]
    f
  })
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}


test_that("it reproduces the published one-atmosphere check values", {
  # UNESCO (1983), density at p = 0.
  expect_equal(eos80_density(5, 0), 999.96675, tolerance = 1e-6)
  expect_equal(eos80_density(5, 35), 1027.67547, tolerance = 1e-6)
  expect_equal(eos80_density(25, 0), 997.04796, tolerance = 1e-6)
  expect_equal(eos80_density(25, 35), 1023.34306, tolerance = 1e-6)
})

test_that("fresh water is densest near 4 degrees", {
  # The anomaly that makes lakes freeze from the top down. A polynomial with a
  # mistyped coefficient will usually lose it.
  t <- seq(0, 12, by = 0.05)
  peak <- t[which.max(eos80_density(t, 0))]
  expect_equal(peak, 3.98, tolerance = 0.05)
})

test_that("seawater has no density maximum in the ocean range", {
  # At 35 PSU the maximum has been pushed below the freezing point, so density
  # falls monotonically with temperature.
  t <- seq(0, 30, by = 0.5)
  expect_true(all(diff(eos80_density(t, 35)) < 0))
})

test_that("saltier is denser at fixed temperature", {
  s <- seq(0, 40, by = 0.5)
  expect_true(all(diff(eos80_density(10, s)) > 0))
})

test_that("sigma-theta is density minus 1000", {
  env <- ts_field(t = c(5, 25), s = c(35, 35))

  sig <- potential_density(env)
  abs_rho <- potential_density(env, sigma = FALSE)

  expect_equal(sig$sigma_theta, abs_rho$density - 1000)
  expect_equal(unique(round(sig$sigma_theta, 4)),
               round(c(1027.67547, 1023.34306) - 1000, 4))
})

test_that("the column is named for what it holds", {
  env <- ts_field(t = c(5, 25), s = c(35, 35))

  expect_true("sigma_theta" %in% names(potential_density(env)))
  expect_true("density" %in% names(potential_density(env, sigma = FALSE)))
  expect_true("rho" %in% names(potential_density(env, name = "rho")))
})

test_that("cold and fresh can be lighter than warm and salty", {
  # The reason temperature alone is a poor stand-in for density: the Scotian
  # Shelf inflow case, cold and fresh arriving over warmer saltier water.
  inflow <- eos80_density(4, 31)
  resident <- eos80_density(10, 35)
  expect_lt(inflow, resident)
})

test_that("salinity given as a fraction is caught rather than silently wrong", {
  # 0.035 instead of 35 is the classic unit error, and produces a number close
  # to fresh water that looks entirely plausible.
  env <- ts_field(t = c(10, 10), s = c(0.035, 0.035))

  # It sits inside the fitted range, so only a units check can catch it.
  expect_warning(potential_density(env), "PSU")
})

test_that("genuinely fresh water still computes, having been flagged", {
  env <- ts_field(t = c(10, 10), s = c(0.035, 0.035))

  suppressWarnings(result <- potential_density(env))

  # Near fresh water, which is the point: the number is right for the input
  # given, and the warning is about whether the input was meant.
  expect_equal(unique(round(result$sigma_theta, 1)), -0.3)
})

test_that("an impossible temperature warns", {
  env <- ts_field(t = c(10, 300), s = c(35, 35))
  expect_warning(potential_density(env), "temperature")
})

test_that("values inside the fitted range do not warn", {
  env <- ts_field(t = c(2, 18), s = c(31, 35))
  expect_no_warning(potential_density(env))
})

test_that("missing values propagate rather than erroring", {
  env <- ts_field(t = c(5, 25), s = c(35, 35))
  env$SST[1] <- NA

  result <- potential_density(env)

  expect_true(is.na(result$sigma_theta[1]))
  expect_false(any(is.na(result$sigma_theta[-1])))
})

test_that("other temperature and salinity columns can be named", {
  env <- ts_field(t = c(5, 25), s = c(35, 35))
  env$BOTT <- env$SST - 2
  env$BOTS <- env$SSS + 0.5

  result <- potential_density(env, temperature = "BOTT", salinity = "BOTS",
                              name = "sigma_bottom")

  expect_equal(result$sigma_bottom, eos80_density(env$BOTT, env$BOTS) - 1000)
})

test_that("a missing column is an error naming it", {
  env <- ts_field(t = 5, s = 35)
  env$SSS <- NULL
  expect_error(potential_density(env), "SSS")
})
