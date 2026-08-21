#' How long water stays in a place
#'
#' Releases a particle at every point inside a box and advects it until it
#' leaves, giving the time the water starting there remained. Long residence
#' means a retentive place — a gyre, a basin, the lee of a bank — where anything
#' with a life stage measured in weeks can complete it without being swept out.
#' Short residence means water is passing through.
#'
#' This is the quantity behind a great deal of plankton distribution that
#' concentration alone cannot explain. *Calanus* accumulates in the deep basins
#' of the Gulf of Maine partly because the circulation holds water there long
#' enough, and a covariate that says so is closer to the mechanism than one that
#' says the water is currently cold.
#'
#' @section Censoring, which is not optional to think about:
#' A particle still inside the box when `max_days` runs out has a residence time
#' of *at least* that, not equal to it. The column records `max_days` for those,
#' which is right-censored data, and treating it as a measurement will bias
#' every summary of it downwards — the more retentive the site, the worse.
#'
#' The function warns with the fraction censored. If that fraction is large the
#' honest options are to raise `max_days`, or to treat the column as censored
#' and model it accordingly, or to use it only as an ordering rather than as a
#' duration. What it is not safe to do is average it.
#'
#' @section Resolution:
#' Membership is tested once per step, so a residence time is the first check
#' after the particle actually left: at or above the true value, by less than
#' `step_hours`. That is the resolution of the answer, and on a short window it
#' is a large share of it: a 6-hour step resolves a 2-day residence to about one
#' part in eight. Shorten the step when the residences being compared are
#' themselves short.
#'
#' @section Cost:
#' One particle per point per time step, integrated up to `max_days` at
#' `step_hours`. That is the same order of work as [ftle()] and rather more of
#' it, since the window is usually longer. Start with a short `max_days` on a
#' small domain.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param box named list with `xmin`, `xmax`, `ymin`, `ymax`, the region water
#'   is considered to be residing in. Defaults to the extent of the data, which
#'   measures retention by the domain itself and is rarely what you want
#' @param u name of the eastward velocity column, in m/s
#' @param v name of the northward velocity column, in m/s
#' @param max_days longest residence measured; particles still inside at this
#'   point are censored here
#' @param step_hours integration step
#' @param direction `"forward"` for how long water will stay, `"backward"` for
#'   how long it has already been there
#' @param name name for the new column
#' @return `env_dat` with a residence time column, in days. `NA` outside the box
#'   and for particles that left the velocity field without leaving the box
#' @examples
#' \dontrun{
#' # Retention in the eastern Gulf of Maine.
#' env <- residence_time(env, box = eastern_gom_box(), max_days = 60)
#' }
#' @seealso [ftle()], [fsle()], [box_anomaly()]
#' @export
residence_time <- function(env_dat, box = NULL, u = "UO", v = "VO",
                           max_days = 60, step_hours = 6,
                           direction = c("forward", "backward"), name = NULL) {
  direction <- match.arg(direction)
  resolve_vars(env_dat, c(u, v))

  if (max_days <= 0) {
    stop("max_days must be positive; use `direction` to look into the past.",
         call. = FALSE)
  }
  if (step_hours <= 0 || step_hours > max_days * 24) {
    stop("step_hours must be positive and no longer than max_days.",
         call. = FALSE)
  }

  box <- if (is.null(box)) extent_box(env_dat) else check_box(box)

  steps <- time_steps(env_dat)
  step_times <- step_time_days(steps)
  velocity <- lapply(seq_len(nrow(steps)), function(i) {
    rasterize_step(env_dat[step_rows(env_dat, steps[i, ]), ], c(u, v))
  })

  sign <- if (direction == "backward") -1 else 1
  step_days <- step_hours / 24
  n_steps <- ceiling(max_days / step_days)

  result <- rep(NA_real_, nrow(env_dat))
  censored <- 0L
  escaped <- 0L
  released <- 0L

  for (i in seq_len(nrow(steps))) {
    rows <- step_rows(env_dat, steps[i, ])
    seeds <- sf::st_coordinates(env_dat[rows, ])
    inside_start <- inside_box(seeds, box)
    if (!any(inside_start)) next

    active <- which(inside_start)
    positions <- seeds[active, , drop = FALSE]
    elapsed <- rep(NA_real_, length(active))
    time <- step_times[i]
    released <- released + length(active)

    for (k in seq_len(n_steps)) {
      # One RK4 step, using the same integrator as ftle() rather than a second
      # copy of it, so the two cannot drift apart.
      positions <- advect(positions, start_time = time,
                          duration = sign * step_days, step_days = step_days,
                          velocity = velocity, times = step_times)
      time <- time + sign * step_days

      still_here <- inside_box(positions, box)
      lost <- !stats::complete.cases(positions)
      # Out of the box is an answer; out of the velocity field is not.
      done <- which(!still_here & !lost & is.na(elapsed))
      elapsed[done] <- k * step_days
      gone <- which(lost & is.na(elapsed))
      elapsed[gone] <- -Inf

      if (all(!is.na(elapsed))) break
    }

    escaped <- escaped + sum(is.infinite(elapsed))
    still_inside <- is.na(elapsed)
    censored <- censored + sum(still_inside)
    elapsed[still_inside] <- max_days
    elapsed[is.infinite(elapsed)] <- NA_real_

    result[rows[active]] <- elapsed
  }

  warn_residence(released, censored, escaped, max_days)
  env_dat[[name %||% paste0(direction, "_residence")]] <- result
  env_dat
}

#' Whether each position lies inside a box
#'
#' @param positions two-column matrix of lon/lat
#' @param box a validated box
#' @return logical vector; `FALSE` where the position is `NA`
#' @keywords internal
inside_box <- function(positions, box) {
  ok <- stats::complete.cases(positions)
  ok & positions[, 1] >= box$xmin & positions[, 1] <= box$xmax &
    positions[, 2] >= box$ymin & positions[, 2] <= box$ymax
}

#' The bounding box of the data
#'
#' @param env_dat an `sf` POINT object
#' @return a box list
#' @keywords internal
extent_box <- function(env_dat) {
  bbox <- sf::st_bbox(env_dat)
  list(xmin = unname(bbox["xmin"]), xmax = unname(bbox["xmax"]),
       ymin = unname(bbox["ymin"]), ymax = unname(bbox["ymax"]))
}

#' Report how much of the result is censored or unknown
#'
#' @param released particles released
#' @param censored particles still inside when the window ended
#' @param escaped particles that left the velocity field without leaving the box
#' @param max_days the window
#' @return invisible `NULL`
#' @keywords internal
warn_residence <- function(released, censored, escaped, max_days) {
  if (released == 0) {
    warning("No points fall inside the box, so nothing was released.",
            call. = FALSE)
    return(invisible(NULL))
  }

  fraction <- censored / released
  if (fraction > 0.1) {
    warning(
      format(round(100 * fraction), nsmall = 0), "% of particles (", censored,
      " of ", released, ") were still inside the box after ", max_days,
      " days.",
      "\n  Their residence time is recorded as ", max_days, " but is really ",
      "at least that, so the column is\n  right-censored and its mean is ",
      "biased downwards -- most severely at the most retentive sites,",
      "\n  which is usually the comparison being made. Raise max_days, or ",
      "treat the column as censored\n  rather than as a duration.",
      call. = FALSE
    )
  }

  if (escaped > 0) {
    warning(
      escaped, " of ", released, " particles left the velocity field while ",
      "still inside the box, and are NA.",
      "\n  That happens where the box reaches the edge of the fetched domain: ",
      "the particle went\n  somewhere the data does not cover, so whether it ",
      "was still residing is unknowable.",
      "\n  Fetch a larger bounding box, or move the box inside the domain.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
