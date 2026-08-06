# A domain too small to hold a trajectory is the commonest reason FTLE and FSLE
# come back empty. Every individual NA is correct, so what is tested here is that
# the emptiness is explained rather than silently returned.

test_that("FTLE warns, with numbers, when every particle leaves the domain", {
  # 0.2 m/s over 14 days is ~240 km. The box is about 80 km across, so nothing
  # survives the window.
  env <- make_flow(function(lon, lat) cbind(rep(-0.2, length(lon)), rep(0, length(lon))),
                   lon = seq(-70, -69, by = 0.05), lat = seq(42, 43, by = 0.05))

  expect_warning(result <- ftle(env, integration_days = 14),
                 "returned no values at all")
  expect_true(all(is.na(result$backward_ftle)))

  # The message has to carry the numbers that make it actionable: how far the
  # flow reaches, how big the box is, and how many particles were lost.
  message <- tryCatch(ftle(env, integration_days = 14), warning = conditionMessage)
  expect_match(message, "14-day integration")
  expect_match(message, "median speed")
  expect_match(message, "particles left the velocity field")
  expect_match(message, "upstream edge")
})

test_that("the named edge follows the direction of integration", {
  env <- make_flow(function(lon, lat) cbind(rep(-0.2, length(lon)), rep(0, length(lon))),
                   lon = seq(-70, -69, by = 0.05), lat = seq(42, 43, by = 0.05))

  backward <- tryCatch(ftle(env, integration_days = 14), warning = conditionMessage)
  forward <- tryCatch(ftle(env, integration_days = 14, direction = "forward"),
                      warning = conditionMessage)

  expect_match(backward, "upstream edge")
  expect_match(forward, "downstream edge")
})

test_that("a domain that comfortably holds the trajectory is silent", {
  # The warning must not fire on ordinary use. Losing a margin to the edge is
  # normal; losing essentially everything is not.
  env <- make_flow(function(lon, lat) cbind(rep(0.02, length(lon)), rep(0, length(lon))),
                   lon = seq(-70, -68, by = 0.05), lat = seq(42, 43, by = 0.05))

  expect_no_warning(result <- ftle(env, integration_days = 2, step_hours = 12))
  expect_gt(sum(!is.na(result$backward_ftle)), 0)
})

test_that("FSLE blames the domain when parcels are lost before separating", {
  # Fast uniform flow: pairs never separate, but they leave the box first, so
  # the binding constraint is max_days against the size of the box.
  env <- make_flow(function(lon, lat) cbind(rep(-0.2, length(lon)), rep(0, length(lon))),
                   lon = seq(-70, -69, by = 0.05), lat = seq(42, 43, by = 0.05))

  message <- tryCatch(fsle(env, final_separation = 50, max_days = 60),
                      warning = conditionMessage)

  expect_match(message, "drifted out of the velocity field")
  expect_match(message, "binding constraint is max_days")
})

test_that("FSLE blames the strain when parcels stay put and never separate", {
  # Almost motionless: nothing leaves the box, and nothing separates either.
  env <- make_flow(function(lon, lat) cbind(rep(1e-4, length(lon)), rep(0, length(lon))),
                   lon = seq(-70, -68, by = 0.05), lat = seq(42, 43, by = 0.05))

  message <- tryCatch(fsle(env, final_separation = 50, max_days = 5),
                      warning = conditionMessage)

  expect_match(message, "never reached the 50 km target")
  expect_match(message, "never separated")
})

test_that("the two FSLE diagnoses are distinguished, not merged", {
  lost_case <- make_flow(function(lon, lat) cbind(rep(-0.2, length(lon)), rep(0, length(lon))),
                         lon = seq(-70, -69, by = 0.05), lat = seq(42, 43, by = 0.05))
  still_case <- make_flow(function(lon, lat) cbind(rep(1e-4, length(lon)), rep(0, length(lon))),
                          lon = seq(-70, -68, by = 0.05), lat = seq(42, 43, by = 0.05))

  lost <- tryCatch(fsle(lost_case, max_days = 60), warning = conditionMessage)
  still <- tryCatch(fsle(still_case, max_days = 5), warning = conditionMessage)

  # Each should give the advice that fits its own cause, and not the other's.
  expect_false(grepl("never separated", lost))
  expect_false(grepl("binding constraint is max_days", still))
})

test_that("a partially empty field is not warned about", {
  # Warning whenever any margin is lost would fire on nearly every real call and
  # train the user to ignore it.
  expect_false(warned(warn_lagrangian("ftle", 0.5, "d", "a")))
  expect_false(warned(warn_lagrangian("ftle", 0.89, "d", "a")))
  expect_true(warned(warn_lagrangian("ftle", 0.9, "d", "a")))
  expect_true(warned(warn_lagrangian("ftle", 1, "d", "a")))
})

test_that("the headline distinguishes empty from nearly empty", {
  expect_match(tryCatch(warn_lagrangian("ftle", 1, "d", "a"), warning = conditionMessage),
               "no values at all")
  expect_match(tryCatch(warn_lagrangian("ftle", 0.95, "d", "a"), warning = conditionMessage),
               "almost nothing")
})
