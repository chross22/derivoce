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
#' The order is deliberate: the trend is removed first, then the seasonal cycle
#' is taken from the detrended values. Doing it the other way lets a trend
#' sitting unevenly across the calendar leak into the seasonal term — with a
#' record starting in July and warming throughout, the later-summer months carry
#' more warm years than the winter ones and the "seasonal cycle" acquires a step
#' that is really the trend.
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

  warn_decompose(short, thin, length(cell_rows), degree)
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

  centre <- mean(y[usable])
  fit <- stats::lm(y ~ stats::poly(t, degree, raw = TRUE),
                   data = data.frame(y = y, t = t), na.action = stats::na.exclude)
  fitted <- stats::predict(fit, newdata = data.frame(t = t))

  trend <- fitted - mean(fitted[usable])
  detrended <- y - centre - trend

  # The seasonal term is only meaningful once a month has more than one year
  # behind it; otherwise it is that month's own value and the residual is zero
  # by construction.
  per_month <- table(month[usable])
  if (length(per_month) == 0 || max(per_month) < 2) out$too_thin <- 1L

  monthly <- tapply(detrended[usable], month[usable], mean)
  seasonal <- as.numeric(monthly[match(month, names(monthly))])
  seasonal <- seasonal - mean(seasonal[usable], na.rm = TRUE)

  out$trend <- trend
  out$seasonal <- seasonal
  out$residual <- detrended - seasonal
  if (degree == 1) {
    # Slope per year rather than per day, which is the unit anyone reads.
    out$slope <- rep(unname(stats::coef(fit)[2]) * 365.25, length(y))
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
      thin, " of ", total, " cell(s) have at most one value per calendar ",
      "month, so the seasonal term is\n  that month's own value and the ",
      "residual is zero by construction.",
      "\n  A seasonal cycle needs several years. Use components = \"trend\" ",
      "on a record this short.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
