#' Rolling summary of a covariate over a trailing window
#'
#' Summarises each location's recent history: the mean of the last three months,
#' the variability of the last year, the coldest step in the last six. Where
#' [integrate_covariate()] accumulates a total, this describes the distribution
#' the total came from, which is often the more useful covariate — a mean and a
#' standard deviation say different things about a place than their sum does.
#'
#' The window is **trailing and inclusive**: it ends at the current step and
#' includes it. A three-month mean at March covers January, February and March.
#'
#' @section Steps or calendar time:
#' `by = "step"` counts positions in the record and `by = "day"`, `"month"` or
#' `"year"` count calendar time, exactly as in [lag_covariate()]. The two agree
#' until the record has a gap and then disagree silently: on a monthly series
#' missing April, a three-*step* window at June covers March, May and June,
#' while a three-*month* window covers April, May and June and finds only two of
#' them.
#'
#' Which is right depends on the question. "The mean of the last three months"
#' is a statement about the ocean and wants `by = "month"`. "The mean of the
#' last three observations" is a statement about the record.
#'
#' @section Windows that are not full:
#' Early steps have less history behind them than the window asks for, and a
#' location can be absent from some of the steps in it. `min_obs` sets how many
#' values a window must actually contain before it is summarised; below that the
#' result is `NA` rather than a mean of whatever happened to be there. The
#' default of 1 is permissive, so the first steps of a record get a summary of a
#' short window rather than nothing. Raise it if a partial window would mislead.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param vars covariate columns, or `NULL` for all numeric ones
#' @param n length of the window, in `by` units, including the current step
#' @param by `"step"` to count positions in the record, or `"day"`, `"month"`,
#'   `"year"` to count calendar time
#' @param stat one or more of `"mean"`, `"sd"`, `"min"`, `"max"`, `"sum"`,
#'   `"median"`, `"range"`
#' @param min_obs fewest non-missing values a window must hold to be summarised
#' @param suffix one per `stat`, or `NULL` to name them automatically
#' @return `env_dat` with one column per covariate per statistic
#' @examples
#' \dontrun{
#' # Conditions over the season leading up to each observation.
#' env <- rolling_covariate(env, "SST", n = 3, by = "month")
#'
#' # How variable it has been, which is a different covariate from how warm.
#' env <- rolling_covariate(env, "SST", n = 12, by = "month", stat = "sd")
#' }
#' @seealso [integrate_covariate()], [lag_covariate()], [cell_anomaly()]
#' @export
rolling_covariate <- function(env_dat, vars = NULL, n = 3,
                              by = c("step", "day", "month", "year"),
                              stat = c("mean", "sd", "min", "max", "sum",
                                       "median", "range"),
                              min_obs = 1L, suffix = NULL) {
  by <- match.arg(by)
  stat <- match.arg(stat, several.ok = TRUE)
  vars <- resolve_vars(env_dat, vars, kind = "temporal")

  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1 || n != round(n)) {
    stop("`n` must be a single whole number of at least 1.", call. = FALSE)
  }
  if (!is.null(suffix) && length(suffix) != length(stat)) {
    stop("`suffix` must have one entry per statistic, or be NULL to name them ",
         "automatically.", call. = FALSE)
  }
  unit <- if (by == "step") "" else by
  suffixes <- suffix %||% paste0("_", stat, n, unit)

  steps <- time_steps(env_dat)
  location <- location_key(env_dat)
  step_index <- match_step_index(env_dat, steps)

  for (v in vars) {
    results <- lapply(stat, function(s) rep(NA_real_, nrow(env_dat)))
    names(results) <- stat

    for (i in seq_len(nrow(steps))) {
      current <- which(step_index == i)
      if (length(current) == 0) next
      contributing <- window_steps(steps, i, n, by)

      window <- matrix(NA_real_, nrow = length(current),
                       ncol = length(contributing))
      for (k in seq_along(contributing)) {
        earlier <- which(step_index == contributing[k])
        window[, k] <- env_dat[[v]][earlier][match(location[current],
                                                   location[earlier])]
      }

      counted <- rowSums(!is.na(window))
      for (s in stat) {
        summarised <- apply(window, 1, function(z) summarise_window(z, s))
        summarised[counted < min_obs] <- NA_real_
        results[[s]][current] <- summarised
      }
    }

    for (k in seq_along(stat)) {
      env_dat[[paste0(v, suffixes[k])]] <- results[[stat[k]]]
    }
  }
  env_dat
}

#' Apply one summary to a window, ignoring missing values
#'
#' `min()` and `max()` of an all-missing window would be `Inf` with a warning
#' rather than `NA`, and `sd()` of a single value is `NA` already, so the empty
#' case is handled once here.
#'
#' @param z the values in the window
#' @param stat the statistic to apply
#' @return a single number
#' @keywords internal
summarise_window <- function(z, stat) {
  z <- z[!is.na(z)]
  if (length(z) == 0) return(NA_real_)
  switch(
    stat,
    mean = mean(z),
    sd = if (length(z) < 2) NA_real_ else stats::sd(z),
    min = min(z),
    max = max(z),
    sum = sum(z),
    median = stats::median(z),
    range = max(z) - min(z)
  )
}

#' Which steps fall in the trailing window ending at a given step
#'
#' The calendar cases are done on the same month counter as [lag_source_step()],
#' so that a window and a lag of the same size agree about what a month is.
#'
#' @param steps a time-step table from [time_steps()]
#' @param i the step the window ends at
#' @param n window length in `by` units, including step `i`
#' @param by `"step"`, `"day"`, `"month"`, or `"year"`
#' @return the indices of the contributing steps
#' @keywords internal
window_steps <- function(steps, i, n, by) {
  if (identical(by, "step")) return(seq(max(1L, i - n + 1L), i))

  if (identical(by, "day")) {
    # In days with an hour fraction where there is one, so on a sub-daily
    # record the window trails from this instant rather than sweeping in
    # later hours of the current day.
    when <- step_time_days(steps)
    return(which(when <= when[i] & when > when[i] - n))
  }

  months_back <- if (identical(by, "year")) n * 12 else n
  counter <- steps$YEAR * 12 + (steps$MONTH - 1)
  which(counter <= counter[i] & counter > counter[i] - months_back)
}
