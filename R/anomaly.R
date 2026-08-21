#' Anomaly of a covariate against each cell's own history
#'
#' Subtracts, from every value, the mean of that same grid cell. The result is a
#' **map per time step** saying where conditions departed from what that place
#' is usually like, rather than one number per step.
#'
#' This is the cell-wise counterpart of [box_anomaly()]. Where that averages a
#' region into an index, this leaves the spatial pattern intact and removes the
#' spatial pattern of the *mean* instead. It is often what makes a covariate
#' usable across a domain with a strong background gradient: 8 degrees is cold
#' for the southern Gulf of Maine and warm for the Scotian Shelf, and a model
#' given raw temperature has to learn that geography before it can use the
#' departure. An anomaly hands it over directly.
#'
#' @section Standardising:
#' `standardize = TRUE` divides by the cell's standard deviation as well,
#' giving a z-score. That makes departures comparable between places with
#' different variability — a 1 degree anomaly is unremarkable on the shelf and
#' extreme in the deep basin — which is what you want when a single coefficient
#' has to apply across the whole domain.
#'
#' It also throws away the magnitude. If the question is how much warmer, not
#' how unusual, leave it off.
#'
#' @section Detrending:
#' `detrend = TRUE` removes a fitted linear trend as well. This matters in a
#' warming shelf sea: an anomaly that still contains the trend largely encodes
#' *which year it is*, and a model given it will fit the trend and appear to
#' have learned something about temperature. What is left after detrending says
#' whether conditions were unusual **for their year**, which is usually the
#' ecological question.
#'
#' The trend and the seasonal cycle are estimated in one fit rather than
#' sequentially, for the reasons set out in [decompose_covariate()] — removed
#' one after the other, each absorbs part of the other. With
#' `reference = "climatology"` both are taken out and the result is that
#' function's residual; with `reference = "record"` only the mean and the trend
#' go, and the seasonal cycle stays in.
#'
#' Detrending is not free of assumptions. A linear trend fitted to a short
#' record can absorb genuine low-frequency variability — a decade of a
#' multidecadal oscillation looks like a trend — so the residual is a departure
#' from a fitted line, not from a known baseline.
#'
#' @section Why this needs several years:
#' With `reference = "climatology"` the mean is taken over the values sharing a
#' calendar month, so a cell's January is compared with its other Januaries.
#' Given a single year of data there is exactly one January per cell, the mean
#' is that value, and **every anomaly is identically zero**. That is a silent
#' failure — a column of zeroes is a perfectly plausible-looking covariate — so
#' it warns.
#'
#' The same applies to `standardize`, which needs at least two values per group
#' to have a standard deviation at all.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param vars covariate columns, or `NULL` for all numeric ones
#' @param reference `"climatology"` (the default) removes a separate mean per
#'   cell per calendar month, leaving departures from the usual conditions for
#'   the time of year. `"record"` removes one mean per cell over the whole
#'   series, which leaves the seasonal cycle in
#' @param standardize divide by the standard deviation of the same group, giving
#'   a z-score
#' @param detrend also remove a fitted linear trend from each cell, so what is
#'   left is the departure after both the seasonal cycle and the long-term
#'   change are gone
#' @param suffix appended to each covariate name to make the new column. Defaults
#'   to `"_anom"`, or `"_z"` when standardising
#' @return `env_dat` with one anomaly column per covariate, in the units of the
#'   covariate, or dimensionless when standardised
#' @examples
#' \dontrun{
#' # Departures from the usual conditions for the month, cell by cell.
#' env <- cell_anomaly(env, "SST")
#'
#' # Comparable across a domain with very different variability.
#' env <- cell_anomaly(env, "SST", standardize = TRUE)
#' }
#' @seealso [box_anomaly()], [marine_heatwave()], [lag_covariate()]
#' @export
cell_anomaly <- function(env_dat, vars = NULL,
                         reference = c("climatology", "record"),
                         standardize = FALSE, detrend = FALSE, suffix = NULL) {
  reference <- match.arg(reference)
  vars <- resolve_vars(env_dat, vars, kind = "temporal")
  suffix <- suffix %||% if (standardize) "_z" else "_anom"

  groups <- anomaly_groups(env_dat, reference)
  if (!detrend) warn_thin_climatology(groups, reference, standardize)

  for (v in vars) {
    values <- env_dat[[v]]
    centred <- if (detrend) {
      detrended_anomaly(env_dat, values, reference)
    } else {
      values - group_stat(values, groups, mean)
    }
    # Within a group the anomaly differs from the values by a constant, so for
    # the undetrended case this is the same divisor as the values would give.
    if (standardize) centred <- centred / group_stat(centred, groups, stats::sd)
    env_dat[[paste0(v, suffix)]] <- centred
  }
  env_dat
}

#' Anomaly with the long-term trend removed as well
#'
#' Uses the same joint fit as [decompose_covariate()] rather than removing a
#' trend and a climatology one after the other, because in sequence each absorbs
#' part of the other.
#'
#' @param env_dat an `sf` POINT object
#' @param values the covariate
#' @param reference `"climatology"` to drop the seasonal cycle too, `"record"`
#'   to keep it
#' @return a numeric vector the same length as `values`
#' @keywords internal
detrended_anomaly <- function(env_dat, values, reference) {
  cells <- location_key(env_dat)
  elapsed <- elapsed_days(env_dat)
  months <- env_dat$MONTH

  out <- rep(NA_real_, length(values))
  thin <- 0L
  cell_rows <- split(seq_along(values), cells)

  for (rows in cell_rows) {
    piece <- decompose_one(values[rows], elapsed[rows], months[rows], degree = 1)
    thin <- thin + piece$too_thin
    out[rows] <- if (reference == "climatology") {
      piece$residual
    } else {
      # Keep the seasonal cycle, drop only the mean and the trend.
      piece$seasonal + piece$residual
    }
  }

  if (reference == "climatology" && thin > 0) {
    warning(
      thin, " of ", length(cell_rows), " cell(s) cannot support a seasonal ",
      "term, which needs at least two\n  calendar months present and one of ",
      "them seen in more than one year. For those the\n  anomaly is detrended ",
      "but not deseasonalised, so the seasonal cycle is still in it.",
      call. = FALSE
    )
  }
  out
}

#' Grouping used to form an anomaly
#'
#' @param env_dat an `sf` POINT object
#' @param reference `"climatology"` or `"record"`
#' @return a character vector, one group label per row
#' @keywords internal
anomaly_groups <- function(env_dat, reference) {
  cells <- location_key(env_dat)
  if (reference == "record") return(cells)
  paste(cells, env_dat$MONTH, sep = "|")
}

#' Apply a summary within each group and broadcast it back
#'
#' `stats::ave()` would drop rows where the covariate is `NA`, so the summary is
#' taken explicitly with `na.rm` and matched back by group.
#'
#' @param values a numeric vector
#' @param groups a grouping vector of the same length
#' @param fun a summary function
#' @return a numeric vector the same length as `values`
#' @keywords internal
group_stat <- function(values, groups, fun) {
  summary <- tapply(values, groups, function(z) {
    z <- z[!is.na(z)]
    if (length(z) == 0) NA_real_ else fun(z)
  })
  as.numeric(summary[match(groups, names(summary))])
}

#' Warn when there is not enough history to form the requested anomaly
#'
#' @param groups the grouping vector
#' @param reference `"climatology"` or `"record"`
#' @param standardize whether a standard deviation is also needed
#' @return invisible `NULL`, called for the warning
#' @keywords internal
warn_thin_climatology <- function(groups, reference, standardize) {
  per_group <- table(groups)
  most <- if (length(per_group) == 0) 0 else max(per_group)

  if (most < 2) {
    warning(
      "Every ", if (reference == "climatology") "cell-month" else "cell",
      " has a single value, so its mean is that value and every anomaly is ",
      "exactly zero.",
      if (reference == "climatology")
        "\n  A monthly climatology needs several years: with one year there is "
      else "\n  ",
      if (reference == "climatology")
        "one January per cell to average."
      else "one time step per cell to average.",
      "\n  Fetch a longer series, or use reference = \"record\" to compare ",
      "each cell against\n  its own whole-series mean instead.",
      call. = FALSE
    )
    return(invisible(NULL))
  }

  if (standardize && min(per_group) < 2) {
    thin <- sum(per_group < 2)
    warning(
      "standardize = TRUE needs at least two values per group, and ", thin,
      " of ", length(per_group), " have one.",
      "\n  Those become NA, because a single value has no standard deviation.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
