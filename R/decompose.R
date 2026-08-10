#' Split a covariate into trend, seasonal cycle, and what is left
#'
#' Separates each cell's series into the parts it is made of: a long-term trend,
#' a repeating seasonal cycle, and the residual departure from both. The pieces
#' are additive, so
#' `value = mean + trend + seasonal + residual` to numerical precision, and any
#' of them can be used on its own.
#'
#' The separation is worth making because the parts answer different questions
#' and are easy to confuse. In a warming shelf sea, a raw temperature anomaly
#' that still contains the trend largely encodes *which year it is*: a model
#' given it will fit the trend and appear to have learned something about
#' temperature. The residual is the part that says whether this month was warm
#' *for its year and season*, which is usually the ecological question.
#'
#' @section What each part is:
#' * **trend** is the fitted long-term component, centred on zero, in the units
#'   of the covariate. `degree = 1` is a straight line; raise it where the trend
#'   curves, and be aware that a high degree will start absorbing variability
#'   that is not a trend.
#' * **seasonal** is the average departure for each calendar month, taken after
#'   the trend is removed and centred on zero. It repeats every year by
#'   construction, so it cannot represent a seasonal cycle that is itself
#'   changing.
#' * **residual** is everything else, which is what [cell_anomaly()] with
#'   `detrend = TRUE` returns.
#' * **slope** is the rate of change per year, one number per cell, repeated on
#'   its rows. Only defined for `degree = 1`, where the trend is a line and a
#'   single rate describes it.
#'
#' @section What it needs, and what it does when it does not have it:
#' Fitting a trend of `degree` needs at least `degree + 1` time steps for that
#' cell, and a monthly climatology needs several years before a month's mean is
#' anything but that month's one value. A cell with too little history gets `NA`
#' for the parts that cannot be estimated rather than a fit through two points,
#' and the function warns.
#'
#' The trend and the seasonal cycle are estimated **together**, in one fit,
#' rather than one being removed before the other is measured. Doing it in sequence
#' lets whichever goes first absorb part of the other, in both directions.
#'
#' A trend removed first picks up part of the seasonal cycle, because the trend
#' is fitted against elapsed days and calendar months are of unequal length: a
#' twelve-month cycle sampled on the first of each month is not orthogonal to a
#' straight line in days, and a series containing nothing but a seasonal cycle
#' comes out with a spurious trend worth a few percent of its amplitude. A
#' seasonal cycle removed first picks up part of the trend whenever the record
#' does not contain whole years — with a series starting in July and warming
#' throughout, the later months carry more warm years than the earlier ones, and
#' the "cycle" acquires a step that is really the trend.
#'
#' Fitting both at once removes each conditional on the other, so neither
#' failure occurs. The cost is that the seasonal term needs enough observations
#' to afford a parameter per calendar month; where it cannot, it is reported as
#' zero rather than estimated badly, and the function warns.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`
#' @param vars covariate columns, or `NULL` for all numeric ones
#' @param degree polynomial degree of the trend. 1 is a straight line
#' @param components which of `"trend"`, `"seasonal"`, `"residual"` and
#'   `"slope"` to add
#' @param suffix appended to each covariate name, before the component name
#' @return `env_dat` with one column per covariate per component
#' @examples
#' \dontrun{
#' env <- decompose_covariate(env, "SST")
#' # adds SST_trend, SST_seasonal, SST_residual
#'
#' # How fast is each cell warming, in degrees per year?
#' env <- decompose_covariate(env, "SST", components = "slope")
#' }
#' @seealso [cell_anomaly()], [box_anomaly()], [index_series()]
#' @export
decompose_covariate <- function(env_dat, vars = NULL, degree = 1,
                                components = c("trend", "seasonal", "residual",
                                               "slope"),
                                suffix = "_") {
  components <- match.arg(components, several.ok = TRUE)
  vars <- resolve_vars(env_dat, vars, kind = "temporal")

  if (!is.numeric(degree) || length(degree) != 1 || is.na(degree) ||
      degree < 1 || degree != round(degree)) {
    stop("`degree` must be a single whole number of at least 1.", call. = FALSE)
  }
  if ("slope" %in% components && degree != 1) {
    stop("`slope` is only defined for degree = 1, where the trend is a line ",
         "and one rate describes it.\n  Drop \"slope\" from `components`, or ",
         "use degree = 1.", call. = FALSE)
  }

  cells <- location_key(env_dat)
  elapsed <- elapsed_days(env_dat)
  months <- env_dat$MONTH

  short <- 0L
  thin <- 0L
  cell_rows <- split(seq_len(nrow(env_dat)), cells)

  for (v in vars) {
    parts <- lapply(components, function(k) rep(NA_real_, nrow(env_dat)))
    names(parts) <- components

    for (rows in cell_rows) {
      piece <- decompose_one(env_dat[[v]][rows], elapsed[rows], months[rows],
                             degree)
      short <- short + piece$too_short
      thin <- thin + piece$too_thin
      for (k in components) parts[[k]][rows] <- piece[[k]]
    }

    for (k in components) {
      env_dat[[paste0(v, suffix, k)]] <- parts[[k]]
    }
  }

  # Only complain about the seasonal term if it was asked for.
  warn_decompose(short, if ("seasonal" %in% components) thin else 0L,
                 length(cell_rows), degree)
  env_dat
}

#' Decompose a single cell's series
#'
#' @param y the covariate values for one cell
#' @param t elapsed days for those rows
#' @param month calendar month for those rows
#' @param degree polynomial degree of the trend
#' @return a list of components, plus counters for the warnings
#' @keywords internal
decompose_one <- function(y, t, month, degree) {
  empty <- rep(NA_real_, length(y))
  out <- list(trend = empty, seasonal = empty, residual = empty, slope = empty,
              too_short = 0L, too_thin = 0L)

  usable <- !is.na(y) & !is.na(t)
  if (sum(usable) < degree + 1) {
    out$too_short <- 1L
    return(out)
  }

  per_month <- table(month[usable])
  months_present <- length(per_month)
  # A month factor costs one parameter per month beyond the first, and is only
  # informative once some month has more than one year behind it. Below that it
  # would fit every observation exactly and drive the residual to zero.
  afford_seasonal <- months_present > 1 && max(per_month) > 1 &&
    sum(usable) > degree + months_present
  if (!afford_seasonal) out$too_thin <- 1L

  frame <- data.frame(y = y, t = t, month = factor(month))
  form <- if (afford_seasonal) {
    y ~ stats::poly(t, degree, raw = TRUE) + month
  } else {
    y ~ stats::poly(t, degree, raw = TRUE)
  }
  fit <- stats::lm(form, data = frame, na.action = stats::na.exclude)

  # Trend and seasonal cycle are estimated together rather than one after the
  # other. Fitting the trend first lets the seasonal cycle leak into it: months
  # are of unequal length, so a twelve-month cycle sampled on the first of each
  # month is not orthogonal to a straight line in elapsed days, and a pure
  # seasonal series acquires a spurious trend of a few percent of its own
  # amplitude. Estimated jointly, neither can absorb the other.
  pieces <- stats::predict(fit, newdata = frame, type = "terms")
  labels <- colnames(pieces)
  trend_col <- grep("poly", labels, fixed = TRUE)

  out$trend <- as.numeric(pieces[, trend_col])
  out$seasonal <- if (afford_seasonal) {
    as.numeric(pieces[, setdiff(seq_along(labels), trend_col)])
  } else {
    rep(0, length(y))
  }
  out$residual <- y - as.numeric(stats::predict(fit, newdata = frame))

  if (degree == 1) {
    slope <- unname(stats::coef(fit)[grep("poly", names(stats::coef(fit)),
                                          fixed = TRUE)][1])
    # Per year rather than per day, which is the unit anyone reads.
    out$slope <- rep(slope * 365.25, length(y))
  }
  out
}

#' Days elapsed since the first time step
#'
#' Real dates rather than step positions, so a gap in the record leaves a gap in
#' the trend rather than compressing it.
#'
#' @param env_dat an `sf` POINT object
#' @return numeric vector of days, one per row
#' @keywords internal
elapsed_days <- function(env_dat) {
  dates <- as.Date(paste(env_dat$YEAR, env_dat$MONTH, env_dat$DAY, sep = "-"))
  as.numeric(dates - min(dates, na.rm = TRUE))
}

#' Warn about cells with too little history to decompose
#'
#' @param short cells with fewer steps than the trend needs
#' @param thin cells whose seasonal term rests on one year
#' @param total number of cells
#' @param degree the requested polynomial degree
#' @return invisible `NULL`
#' @keywords internal
warn_decompose <- function(short, thin, total, degree) {
  if (short > 0) {
    warning(
      short, " of ", total, " cell(s) have fewer than ", degree + 1,
      " time steps, which is the fewest a degree-", degree,
      " trend can be fitted to.",
      "\n  Those cells are NA throughout rather than carrying a fit through ",
      "as many points as it has\n  parameters.",
      call. = FALSE
    )
  }
  if (thin > 0) {
    warning(
      thin, " of ", total, " cell(s) cannot support a seasonal term, which ",
      "needs at least two calendar\n  months present and at least one of them ",
      "observed in more than one year.",
      "\n  Their seasonal component is reported as zero, and whatever a ",
      "seasonal cycle would have\n  explained stays in the residual. On an ",
      "annual series, or one covering a single year,\n  ask for ",
      "components = \"trend\" instead.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
