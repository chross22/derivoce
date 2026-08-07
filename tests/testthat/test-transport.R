# Transport across a section, checked against flows whose answer is exact.

uniform_flow <- function(east = 1, north = 0) {
  make_flow(function(lon, lat) cbind(rep(east, length(lon)), rep(north, length(lon))),
            lon = seq(-70, -68, by = 0.05), lat = seq(42, 43, by = 0.05),
            months = 1)
}

test_that("transport of a uniform flow equals speed times section length", {
  # 1 m/s eastward through a meridional section of known length. Travelling
  # south to north, the right-hand normal points east, so the flow crosses it
  # head-on and the answer is exactly the length in metres.
  env <- uniform_flow(east = 1)

  result <- section_transport(env, from = c(-69, 42.2), to = c(-69, 42.8))
  expected <- 0.6 * 111320

  got <- unique(stats::na.omit(result$transport))
  expect_equal(length(got), 1)
  expect_lt(abs(got - expected) / expected, 0.01)
})

test_that("transport scales with speed", {
  slow <- section_transport(uniform_flow(east = 1), c(-69, 42.2), c(-69, 42.8))
  fast <- section_transport(uniform_flow(east = 2), c(-69, 42.2), c(-69, 42.8))

  expect_equal(unique(stats::na.omit(fast$transport)),
               2 * unique(stats::na.omit(slow$transport)))
})

test_that("reversing the endpoints reverses the sign", {
  # The sign convention is the easiest thing here to get wrong and the hardest
  # to notice, since the magnitude stays right either way.
  env <- uniform_flow(east = 1)

  forward <- section_transport(env, from = c(-69, 42.2), to = c(-69, 42.8))
  backward <- section_transport(env, from = c(-69, 42.8), to = c(-69, 42.2))

  expect_equal(unique(stats::na.omit(forward$transport)),
               -unique(stats::na.omit(backward$transport)))
})

test_that("a section parallel to the flow carries none of it", {
  env <- uniform_flow(east = 1)

  result <- section_transport(env, from = c(-69.5, 42.5), to = c(-68.5, 42.5))

  expect_lt(abs(unique(stats::na.omit(result$transport))), 1)
})

test_that("a section outside the domain returns NA and says so", {
  env <- uniform_flow()

  expect_warning(result <- section_transport(env, from = c(-10, 42.2),
                                             to = c(-10, 42.8)),
                 "returned no values at all")
  expect_true(all(is.na(result$transport)))
})

test_that("endpoints are validated", {
  env <- uniform_flow()

  expect_error(section_transport(env, from = c(-69, 42), to = c(-69, 42)),
               "no length")
  expect_error(section_transport(env, from = c(-69, 142), to = c(-69, 143)),
               "longitude, latitude")
  expect_error(section_transport(env, from = c(-69, 42), to = "x"), "numeric")
})

test_that("a transposed but legal coordinate warns rather than errors", {
  # c(42, -69) reads longitude 42, latitude -69, which is a real place in the
  # Southern Ocean. It cannot be rejected as invalid, so the only honest
  # response is an empty result with a message that names this as a likely
  # cause. Erroring here would mean refusing legitimate southern-hemisphere work.
  env <- uniform_flow()

  expect_warning(result <- section_transport(env, from = c(42, -69),
                                             to = c(43, -69)),
                 "longitude, latitude")
  expect_true(all(is.na(result$transport)))
})

test_that("the named sections are fixed, documented, and oriented into the Gulf", {
  scotian <- scotian_shelf_inflow_section()
  channel <- northeast_channel_section()

  for (section in list(scotian, channel)) {
    expect_named(section, c("from", "to"))
    expect_length(section$from, 2)
    expect_length(section$to, 2)
  }

  # Both are oriented so that the right-hand normal points into the Gulf of
  # Maine, which is to say northwest. If either endpoint pair is ever edited,
  # this is the check that catches a silent sign flip.
  for (section in list(scotian, channel)) {
    dx <- section$to[1] - section$from[1]
    dy <- section$to[2] - section$from[2]
    normal <- c(dy, -dx)
    expect_lt(normal[1], 0)   # westward component
    expect_gt(normal[2], 0)   # northward component
  }
})

test_that("the named indices run and label their own columns", {
  # Placed over the real sections so the flow is actually sampled.
  env <- make_flow(function(lon, lat) cbind(rep(-0.1, length(lon)), rep(0.1, length(lon))),
                   lon = seq(-67, -65, by = 0.05), lat = seq(42, 44, by = 0.05),
                   months = 1)

  expect_true("scotian_inflow" %in% names(scotian_shelf_inflow(env)))
  expect_true("channel_inflow" %in% names(northeast_channel_inflow(env)))
})

test_that("transport is one value per time step, broadcast to every row", {
  # It describes the section, not the cell, so a spatial derivative of it is
  # meaningless and horizontal_gradient() should say so.
  env <- make_flow(function(lon, lat) cbind(rep(1, length(lon)), rep(0, length(lon))),
                   lon = seq(-70, -68, by = 0.1), lat = seq(42, 43, by = 0.1),
                   months = 1:2)

  result <- section_transport(env, from = c(-69, 42.2), to = c(-69, 42.8))

  per_step <- tapply(result$transport, result$MONTH, function(z) length(unique(z)))
  expect_true(all(per_step == 1))
  expect_warning(horizontal_gradient(result, "transport"), "Spatially uniform")
})
