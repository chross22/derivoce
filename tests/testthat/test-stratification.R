# Stratification and baroclinic instability.
#
# Both are analytic given the inputs, so these check against hand-computed
# values rather than against themselves. The sign conventions carry the risk:
# N^2 must be positive when denser water lies below, and the depth ordering
# must be enforced rather than quietly sorted, because swapping the levels
# flips the sign and a stably stratified column would be reported as unstable.

strat_field <- function(rho_shallow, rho_deep, u = 0.4, v = 0, u_deep = 0,
                        v_deep = 0, lat = 43, months = 1:2) {
  grid <- expand.grid(x = c(-70, -69.5), y = lat)
  frames <- lapply(months, function(m) {
    f <- grid
    f$YEAR <- 2020L
    f$MONTH <- as.integer(m)
    f$DAY <- 1L
    f$rho_s <- rho_shallow
    f$rho_d <- rho_deep
    f$UO <- u
    f$VO <- v
    f$UO_deep <- u_deep
    f$VO_deep <- v_deep
    f
  })
  sf::st_as_sf(do.call(rbind, frames), coords = c("x", "y"), crs = 4326)
}

# N^2 = (g / rho0) * (rho_deep - rho_shallow) / (z_deep - z_shallow)
expected_n2 <- function(rho_s, rho_d, dz) (9.81 / 1025) * (rho_d - rho_s) / dz


test_that("N2 matches the definition", {
  env <- strat_field(rho_shallow = 24, rho_deep = 26)

  out <- buoyancy_frequency(env, "rho_s", "rho_d", depths = c(0, 100))

  expect_equal(unique(out$N2), expected_n2(24, 26, 100))
})

test_that("denser water below is stable and gives a positive N2", {
  env <- strat_field(rho_shallow = 24, rho_deep = 26)
  expect_true(all(buoyancy_frequency(env, "rho_s", "rho_d",
                                     depths = c(0, 100))$N2 > 0))
})

test_that("denser water above is unstable, reported rather than hidden", {
  env <- strat_field(rho_shallow = 27, rho_deep = 25)

  expect_warning(out <- buoyancy_frequency(env, "rho_s", "rho_d",
                                           depths = c(0, 100)),
                 "denser water above")
  expect_true(all(out$N2 < 0))
})

test_that("stronger stratification over a thinner layer gives a larger N2", {
  weak <- strat_field(rho_shallow = 25.0, rho_deep = 25.2)
  strong <- strat_field(rho_shallow = 24.0, rho_deep = 27.0)

  a <- buoyancy_frequency(weak, "rho_s", "rho_d", depths = c(0, 100))
  b <- buoyancy_frequency(strong, "rho_s", "rho_d", depths = c(0, 100))
  thin <- buoyancy_frequency(strong, "rho_s", "rho_d", depths = c(0, 20))

  expect_gt(unique(b$N2), unique(a$N2))
  expect_gt(unique(thin$N2), unique(b$N2))
})

test_that("only the density difference matters, not its absolute scale", {
  # Sigma-theta and full density differ by 1000 and must give the same answer.
  sigma <- strat_field(rho_shallow = 24, rho_deep = 26)
  full <- strat_field(rho_shallow = 1024, rho_deep = 1026)

  expect_equal(
    buoyancy_frequency(sigma, "rho_s", "rho_d", depths = c(0, 100))$N2,
    buoyancy_frequency(full, "rho_s", "rho_d", depths = c(0, 100))$N2
  )
})

test_that("frequency returns N, and drops the unstable cells", {
  stable <- strat_field(rho_shallow = 24, rho_deep = 26)
  unstable <- strat_field(rho_shallow = 27, rho_deep = 25)

  n <- buoyancy_frequency(stable, "rho_s", "rho_d", depths = c(0, 100),
                          frequency = TRUE)
  expect_equal(unique(n$N), sqrt(expected_n2(24, 26, 100)))

  suppressWarnings(
    bad <- buoyancy_frequency(unstable, "rho_s", "rho_d", depths = c(0, 100),
                              frequency = TRUE)
  )
  expect_true(all(is.na(bad$N)))
})

test_that("the depths must be given shallower first", {
  env <- strat_field(rho_shallow = 24, rho_deep = 26)

  expect_error(buoyancy_frequency(env, "rho_s", "rho_d", depths = c(100, 0)),
               "shallower first")
  expect_error(buoyancy_frequency(env, "rho_s", "rho_d", depths = c(50, 50)),
               "no layer between")
  expect_error(buoyancy_frequency(env, "rho_s", "rho_d", depths = 100),
               "two numbers")
})

test_that("the column can be named", {
  env <- strat_field(rho_shallow = 24, rho_deep = 26)
  out <- buoyancy_frequency(env, "rho_s", "rho_d", depths = c(0, 100),
                            name = "strat")
  expect_true("strat" %in% names(out))
})


test_that("the Eady growth rate matches the formula", {
  env <- strat_field(rho_shallow = 24, rho_deep = 26, u = 0.4, u_deep = 0)
  env <- buoyancy_frequency(env, "rho_s", "rho_d", depths = c(0, 100))

  out <- eady_growth_rate(env, shallow = c("UO", "VO"),
                          deep = c("UO_deep", "VO_deep"), depths = c(0, 100))

  shear <- 0.4 / 100
  f <- 2 * 7.2921e-5 * sin(43 * pi / 180)
  brunt <- sqrt(expected_n2(24, 26, 100))
  expect_equal(unique(out$eady_growth), 0.31 * f * shear / brunt * 86400)
})

test_that("per second and per day differ by 86400", {
  env <- strat_field(rho_shallow = 24, rho_deep = 26, u = 0.4)
  env <- buoyancy_frequency(env, "rho_s", "rho_d", depths = c(0, 100))

  day <- eady_growth_rate(env, deep = c("UO_deep", "VO_deep"),
                          depths = c(0, 100))
  sec <- eady_growth_rate(env, deep = c("UO_deep", "VO_deep"),
                          depths = c(0, 100), per = "second")

  expect_equal(day$eady_growth, sec$eady_growth * 86400)
})

test_that("more shear grows faster and more stratification grows slower", {
  base <- strat_field(rho_shallow = 24, rho_deep = 26, u = 0.2)
  sheared <- strat_field(rho_shallow = 24, rho_deep = 26, u = 0.8)
  stratified <- strat_field(rho_shallow = 20, rho_deep = 30, u = 0.2)

  rate <- function(e) {
    e <- buoyancy_frequency(e, "rho_s", "rho_d", depths = c(0, 100))
    unique(eady_growth_rate(e, deep = c("UO_deep", "VO_deep"),
                            depths = c(0, 100))$eady_growth)
  }

  expect_gt(rate(sheared), rate(base))
  expect_lt(rate(stratified), rate(base))
})

test_that("shear direction does not matter, only its magnitude", {
  east <- strat_field(rho_shallow = 24, rho_deep = 26, u = 0.5, v = 0)
  north <- strat_field(rho_shallow = 24, rho_deep = 26, u = 0, v = 0.5)

  rate <- function(e) {
    e <- buoyancy_frequency(e, "rho_s", "rho_d", depths = c(0, 100))
    unique(eady_growth_rate(e, deep = c("UO_deep", "VO_deep"),
                            depths = c(0, 100))$eady_growth)
  }

  expect_equal(rate(east), rate(north))
})

test_that("no shear means no growth", {
  env <- strat_field(rho_shallow = 24, rho_deep = 26, u = 0.3, u_deep = 0.3)
  env <- buoyancy_frequency(env, "rho_s", "rho_d", depths = c(0, 100))

  out <- eady_growth_rate(env, deep = c("UO_deep", "VO_deep"),
                          depths = c(0, 100))

  expect_equal(unique(out$eady_growth), 0)
})

test_that("an unstable column has no growth rate rather than an enormous one", {
  env <- strat_field(rho_shallow = 27, rho_deep = 25, u = 0.4)
  suppressWarnings(
    env <- buoyancy_frequency(env, "rho_s", "rho_d", depths = c(0, 100))
  )

  expect_warning(out <- eady_growth_rate(env, deep = c("UO_deep", "VO_deep"),
                                         depths = c(0, 100)),
                 "no growth rate")
  expect_true(all(is.na(out$eady_growth)))
})

test_that("the growth rate rises towards the pole through the Coriolis term", {
  rate <- function(lat) {
    e <- strat_field(rho_shallow = 24, rho_deep = 26, u = 0.4, lat = lat)
    e <- buoyancy_frequency(e, "rho_s", "rho_d", depths = c(0, 100))
    unique(eady_growth_rate(e, deep = c("UO_deep", "VO_deep"),
                            depths = c(0, 100))$eady_growth)
  }
  expect_gt(rate(60), rate(30))
})

test_that("velocity columns must be given in pairs", {
  env <- strat_field(rho_shallow = 24, rho_deep = 26)
  env <- buoyancy_frequency(env, "rho_s", "rho_d", depths = c(0, 100))

  expect_error(eady_growth_rate(env, shallow = "UO",
                                deep = c("UO_deep", "VO_deep"),
                                depths = c(0, 100)),
               "two columns")
})

test_that("a missing column is an error naming it", {
  env <- strat_field(rho_shallow = 24, rho_deep = 26)
  expect_error(buoyancy_frequency(env, "rho_s", "NOPE", depths = c(0, 100)),
               "NOPE")
})
