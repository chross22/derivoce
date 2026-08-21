#' How often a place is frontal
#'
#' The fraction of time steps in which a cell's own gradient was sharp enough to
#' count as a front. Where [distance_to_front()] asks how far the nearest front
#' was at one moment, this asks how reliably a front sits in this place at all.
#'
#' The distinction matters because fronts move. A cell that is frontal in one
#' step out of twenty happened to catch a passing filament; a cell that is
#' frontal in fifteen sits on a persistent feature — a shelf-break front, a tidal
#' mixing front, the edge of a plume — and those are the ones that aggregate
#' plankton reliably enough for a predator to learn. An instantaneous distance
#' cannot tell the two apart, and averaging distance over time does not either,
#' because a cell can be near a different transient front every step.
#'
#' @section Choosing the threshold:
#' Inherited from [distance_to_front()], and the choice matters more here. With
#' `scope = "record"` one cutoff applies throughout, so frequency reflects both
#' how often a front is present and whether this part of the domain is gradient-
#' rich at all. With `scope = "step"` each step is cut at its own quantile, so a
#' fixed fraction of cells is frontal in every step and frequency becomes purely
#' about location. The second is usually what "persistence" is meant to mean.
#'
#' @section The window:
#' `n = NULL`, the default, uses the whole record and gives a static map: one
#' value per cell, repeated on every step. That is a description of the domain
#' rather than a covariate that varies with the observation.
#'
#' Giving `n` makes it a trailing window in the manner of [rolling_covariate()],
#' so it varies through the record and can enter a model alongside conditions at
#' the time. Steps where the cell's gradient is undefined — the outermost ring,
#' or a missing value — are left out of the denominator rather than counted as
#' not frontal.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param var the covariate whose gradient defines a front
#' @param threshold gradient magnitude at or above which a cell is frontal. When
#'   `NULL`, taken from `quantile`
#' @param quantile quantile of the gradient used as the threshold
#' @param scope `"record"` for one cutoff across the series, `"step"` for a
#'   cutoff per time step
#' @param per distance unit for the gradient, `"km"` or `"m"`
#' @param n length of the trailing window, or `NULL` for the whole record
#' @param by `"step"`, `"day"`, `"month"` or `"year"`, as in [rolling_covariate()]
#' @param name name for the new column
#' @return `env_dat` with a frequency column between 0 and 1
#' @examples
#' \dontrun{
#' # A static map of where thermal fronts persist.
#' env <- front_frequency(env, "SST", scope = "step")
#'
#' # How frontal this place has been over the preceding year.
#' env <- front_frequency(env, "SST", scope = "step", n = 12, by = "month")
#' }
#' @seealso [distance_to_front()], [rolling_covariate()]
#' @export
front_frequency <- function(env_dat, var, threshold = NULL, quantile = 0.9,
                            scope = c("record", "step"), per = c("km", "m"),
                            n = NULL, by = c("step", "day", "month", "year"),
                            name = NULL) {
  scope <- match.arg(scope)
  per <- match.arg(per)
  by <- match.arg(by)
  resolve_vars(env_dat, var)

  if (!is.null(threshold) && threshold <= 0) {
    stop("threshold must be positive; a gradient magnitude is never negative.",
         call. = FALSE)
  }
  if (is.null(threshold) && (quantile <= 0 || quantile >= 1)) {
    stop("quantile must be strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.null(n) && (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1 ||
                      n != round(n))) {
    stop("`n` must be a single whole number of at least 1, or NULL for the ",
         "whole record.", call. = FALSE)
  }

  frontal <- frontal_indicator(env_dat, var, threshold, quantile, scope, per)

  steps <- time_steps(env_dat)
  location <- location_key(env_dat)
  step_index <- match_step_index(env_dat, steps)

  frequency <- rep(NA_real_, nrow(env_dat))
  for (i in seq_len(nrow(steps))) {
    current <- which(step_index == i)
    if (length(current) == 0) next
    contributing <- if (is.null(n)) seq_len(nrow(steps)) else
      window_steps(steps, i, n, by)

    window <- matrix(NA_real_, nrow = length(current),
                     ncol = length(contributing))
    for (k in seq_along(contributing)) {
      earlier <- which(step_index == contributing[k])
      window[, k] <- frontal[earlier][match(location[current],
                                            location[earlier])]
    }
    # An undefined gradient is not evidence of no front, so it leaves the
    # denominator rather than counting as a zero.
    counted <- rowSums(!is.na(window))
    frequency[current] <- ifelse(counted == 0, NA_real_,
                                 rowSums(window, na.rm = TRUE) / counted)
  }

  env_dat[[name %||% paste0(var, "_front_freq")]] <- frequency
  env_dat
}

#' Whether each row's own cell is frontal
#'
#' Shares its threshold logic with [distance_to_front()], so that a cell counted
#' as frontal here is one that would have had a distance of zero there.
#'
#' @param env_dat an `sf` POINT object
#' @param var the covariate whose gradient defines a front
#' @param threshold explicit gradient cutoff, or `NULL`
#' @param quantile quantile used when `threshold` is `NULL`
#' @param scope `"record"` or `"step"`
#' @param per `"km"` or `"m"`
#' @return numeric vector of 1, 0 and `NA`
#' @keywords internal
frontal_indicator <- function(env_dat, var, threshold, quantile, scope, per) {
  gradient_column <- paste0(var, "__freq_grad")
  with_gradient <- horizontal_gradient(env_dat, var, per = per,
                                       suffix = "__freq_grad")
  gradients <- with_gradient[[gradient_column]]

  record_threshold <- threshold %||%
    stats::quantile(gradients, quantile, na.rm = TRUE, names = FALSE)

  steps <- time_steps(env_dat)
  indicator <- rep(NA_real_, nrow(env_dat))

  for (i in seq_len(nrow(steps))) {
    rows <- step_rows(env_dat, steps[i, ])
    step_gradients <- gradients[rows]

    cutoff <- if (!is.null(threshold) || scope == "record") {
      record_threshold
    } else {
      stats::quantile(step_gradients, quantile, na.rm = TRUE, names = FALSE)
    }
    indicator[rows] <- as.numeric(step_gradients >= cutoff)
  }
  indicator
}
