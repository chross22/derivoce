#' Pull a per-time-step index back out as a series
#'
#' Several functions here compute **one value per time step** and broadcast it
#' onto every row, so the object keeps its shape and can carry on down a pipe:
#' [box_anomaly()], [section_transport()], [eastern_gom_salinity()],
#' [northeast_channel_inflow()] and the rest of the region-scale indices all
#' behave this way. That is right for modelling, where the covariate has to line
#' up with the observations, and wrong for almost everything else. Plotting a
#' 22-year monthly index from a broadcast column means plotting each value a few
#' thousand times; writing one out means exporting a file mostly made of
#' repetition.
#'
#' This collapses them back: one row per time step, in date order, as a plain
#' data frame.
#'
#' @section What it will not do:
#' A column that varies within a time step is a map, not an index, and
#' collapsing it would silently throw away the spatial pattern and keep an
#' arbitrary one of its values. Naming such a column is an error rather than a
#' quiet mean. If a summary of a map per step is what you want, that is a
#' different operation and an explicit one — take the mean yourself, or use
#' [box_anomaly()], which is exactly that with a region attached.
#'
#' With `vars = NULL` the constant-within-step columns are found for you, so
#' passing an object through several index functions and then calling this
#' returns all of them at once.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param vars index columns to extract, or `NULL` to take every column that is
#'   constant within each time step
#' @return a data frame with `YEAR`, `MONTH`, `DAY` and one column per index,
#'   one row per time step, ordered by date. Not an `sf` object: an index has no
#'   location
#' @examples
#' \dontrun{
#' env <- eastern_gom_salinity(env)
#' env <- section_transport(env, from = c(-67.5, 44.5), to = c(-66.0, 43.5))
#'
#' series <- index_series(env)
#' plot(with(series, as.Date(paste(YEAR, MONTH, DAY, sep = "-"))),
#'      series$egom_salinity, type = "l")
#' }
#' @seealso [box_anomaly()], [section_transport()], [derived_indices()]
#' @export
index_series <- function(env_dat, vars = NULL) {
  steps <- time_steps(env_dat)
  step_index <- match_step_index(env_dat, steps)
  available <- covariate_columns(env_dat)

  is_index <- function(v) {
    is.numeric(env_dat[[v]]) && constant_within(env_dat[[v]], step_index)
  }

  if (is.null(vars)) {
    vars <- available[vapply(available, is_index, logical(1))]
    if (length(vars) == 0) {
      stop("No columns are constant within a time step, so there is no index ",
           "to extract.",
           "\n  Region-scale indices such as box_anomaly() and ",
           "section_transport() produce them;\n  a covariate that varies ",
           "across the grid is a map, and collapsing it would discard that.",
           call. = FALSE)
    }
  } else {
    missing <- setdiff(vars, available)
    if (length(missing) > 0) {
      stop("Column(s) not present: ", paste(missing, collapse = ", "),
           "\nAvailable: ", paste(available, collapse = ", "), call. = FALSE)
    }
    varying <- vars[!vapply(vars, is_index, logical(1))]
    if (length(varying) > 0) {
      stop("Column(s) vary within a time step: ",
           paste(varying, collapse = ", "),
           "\n  These are maps rather than indices, and collapsing one to a ",
           "single value per step would\n  keep an arbitrary cell and discard ",
           "the pattern. Summarise it deliberately instead --\n  ",
           "box_anomaly() averages over a named region.", call. = FALSE)
    }
  }

  out <- steps[, intersect(time_columns(), names(steps)),
               drop = FALSE]
  for (v in vars) {
    # Constant within the step, so the first non-missing value represents it.
    out[[v]] <- vapply(seq_len(nrow(steps)), function(i) {
      values <- env_dat[[v]][step_index == i]
      values <- values[!is.na(values)]
      if (length(values) == 0) NA_real_ else values[1]
    }, numeric(1))
  }

  rownames(out) <- NULL
  out
}
