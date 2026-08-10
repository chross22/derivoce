#' Marine heatwaves and cold spells
#'
#' Flags periods when a cell was unusually warm — or unusually cold — for the
#' time of year, and measures how unusual and for how long. Follows the
#' definition of Hobday et al. (2016): a threshold set by a high percentile of
#' the cell's own seasonally varying climatology, an event being a run of
#' consecutive steps above it, and an intensity measured from the climatology
#' rather than from the threshold.
#'
#' The Gulf of Maine is one of the fastest-warming shelf seas on record, so
#' whether an animal was there in an ordinary summer or a heatwave summer is
#' often a more useful covariate than the temperature itself.
#'
#' @section What each measure is:
#' * **event** is `TRUE` where the cell is in an event that met `min_steps`.
#' * **intensity** is the departure from the climatological mean, in the units
#'   of `var`, and is `NA` outside an event. Measured from the climatology, not
#'   from the threshold, so it is comparable between cells whose thresholds
#'   differ.
#' * **category** is 1 to 4 — moderate, strong, severe, extreme — from how many
#'   multiples of the threshold's excess over the climatology the value reaches.
#' * **duration** is the length of the event in time steps, repeated on every
#'   row belonging to it.
#' * **cumulative** is the summed intensity over the whole event, which
#'   distinguishes a long mild event from a short violent one.
#'
#' @section Time steps, not days:
#' Hobday et al. define an event as five or more consecutive **days**. This
#' works in whatever step the data is in, because
#' `datamatch::accessEnvDat()` serves monthly as readily as daily. On daily data
#' `min_steps = 5` and `max_gap = 2` reproduce the published definition. On
#' monthly data they do not — five consecutive months is a far rarer and larger
#' thing — so the defaults here are deliberately permissive and the choice is
#' left to you.
#'
#' Runs are counted over consecutive entries in the sequence of time steps
#' **present in the data**, which has two consequences worth knowing.
#'
#' A cell absent from a step that other cells have breaks that cell's event: it
#' cannot be shown to have stayed warm through a step it has no value for.
#'
#' A step missing from the record *entirely* is a different matter. It is not in
#' the sequence at all, so the steps either side of it are adjacent and an event
#' runs straight through. On a record with holes in it that can join what were
#' two events, and the duration will be in steps you have rather than in elapsed
#' time. Check the series is complete before reading duration as a period.
#'
#' @section Why it needs a long series:
#' A 90th percentile taken over three values is the largest of the three, and
#' every third step will then be a heatwave. The published definition uses a
#' 30-year baseline. This warns when the groups behind the threshold are too
#' thin for the percentile to mean anything, but it cannot tell you that a
#' 10-year baseline is short for your purpose.
#'
#' Because the climatology is computed from the data given, a warming trend
#' within the series raises the threshold and later heatwaves are measured
#' against a warmer baseline. That is a choice, not an oversight: it makes
#' events relative to recent conditions. If you want a fixed baseline, compute
#' the threshold on a subset and apply it.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`
#' @param var the covariate to examine, typically sea surface temperature
#' @param percentile threshold percentile, between 0 and 1. Hobday et al. use
#'   0.9 for heatwaves; the cold-spell equivalent is 0.1, which `direction`
#'   applies for you
#' @param min_steps shortest run of consecutive steps counted as an event. Use
#'   5 on daily data for the published definition
#' @param max_gap runs separated by at most this many non-exceeding steps are
#'   joined into one event. Use 2 on daily data for the published definition
#' @param direction `"warm"` for heatwaves, `"cold"` for cold spells
#' @param measures which of `"event"`, `"intensity"`, `"category"`,
#'   `"duration"` and `"cumulative"` to add
#' @param prefix stem for the new column names. Defaults to `"mhw"`, or `"mcs"`
#'   when `direction = "cold"`
#' @return `env_dat` with one column per requested measure
#' @references
#' Hobday AJ, Alexander LV, Perkins SE, Smale DA, Straub SC, Oliver ECJ,
#' Benthuysen JA, Burrows MT, Donat MG, Feng M, Holbrook NJ, Moore PJ,
#' Scannell HA, Sen Gupta A, Wernberg T (2016). A hierarchical approach to
#' defining marine heatwaves. *Progress in Oceanography* **141**, 227-238.
#' \doi{10.1016/j.pocean.2015.12.014}
#'
#' Hobday AJ, Oliver ECJ, Sen Gupta A, Benthuysen JA, Burrows MT, Donat MG,
#' Holbrook NJ, Moore PJ, Thomsen MS, Wernberg T, Smale DA (2018). Categorizing
#' and naming marine heatwaves. *Oceanography* **31**(2), 162-173.
#' \doi{10.5670/oceanog.2018.205}
#' @examples
#' \dontrun{
#' env <- marine_heatwave(env, "SST")
#'
#' # The published definition, on daily data.
#' env <- marine_heatwave(env, "SST", min_steps = 5, max_gap = 2)
#'
#' # Cold spells instead.
#' env <- marine_heatwave(env, "SST", direction = "cold")
#' }
#' @seealso [cell_anomaly()], [box_anomaly()]
#' @export
marine_heatwave <- function(env_dat, var = "SST", percentile = 0.9,
                            min_steps = 1L, max_gap = 0L,
                            direction = c("warm", "cold"),
                            measures = c("event", "intensity", "category",
                                         "duration", "cumulative"),
                            prefix = NULL) {
  direction <- match.arg(direction)
  measures <- match.arg(measures, several.ok = TRUE)
  resolve_vars(env_dat, var, kind = "temporal")
  prefix <- prefix %||% if (direction == "warm") "mhw" else "mcs"

  if (!is.numeric(percentile) || length(percentile) != 1 ||
      !is.finite(percentile) || percentile <= 0 || percentile >= 1) {
    stop("`percentile` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  # A cold spell is the same computation on the other tail.
  if (direction == "cold" && percentile > 0.5) percentile <- 1 - percentile

  values <- env_dat[[var]]
  groups <- anomaly_groups(env_dat, "climatology")
  warn_thin_threshold(groups, percentile)

  climatology <- group_stat(values, groups, mean)
  threshold <- group_stat(values, groups,
                          function(z) stats::quantile(z, percentile, names = FALSE))

  # Everything below is written warm-side-up; a cold spell is the mirror image,
  # so the sign is flipped once here and the same code serves both.
  sign <- if (direction == "warm") 1 else -1
  excess <- sign * (values - climatology)
  margin <- sign * (threshold - climatology)
  exceeds <- !is.na(excess) & (excess > margin)

  events <- event_runs(env_dat, exceeds, min_steps = min_steps,
                       max_gap = max_gap)

  in_event <- events$id > 0
  intensity <- ifelse(in_event, sign * excess, NA_real_)

  out <- list(
    event = in_event,
    intensity = intensity,
    category = ifelse(in_event, mhw_category(excess, margin), NA_integer_),
    duration = ifelse(in_event, events$duration, NA_integer_),
    cumulative = ifelse(in_event, event_sum(events$id, excess), NA_real_)
  )[measures]

  for (m in names(out)) env_dat[[paste0(prefix, "_", m)]] <- out[[m]]
  env_dat
}

#' Hobday category from how far past the threshold a value reached
#'
#' The categories are multiples of the threshold's own excess over the
#' climatology, so a place with little variability reaches a high category on a
#' smaller absolute departure than a variable one.
#'
#' @param excess departure from the climatology, sign-corrected
#' @param margin the threshold's departure from the climatology, sign-corrected
#' @return an integer vector, 1 to 4
#' @keywords internal
mhw_category <- function(excess, margin) {
  ratio <- excess / margin
  # A zero margin means the threshold sits on the climatology and the ratio is
  # undefined; anything above it is at least moderate.
  ratio[!is.finite(ratio)] <- 1
  pmin(4L, pmax(1L, as.integer(floor(ratio))))
}

#' Label runs of consecutive exceeding time steps, per cell
#'
#' Works in step index rather than row order, so a cell absent from a step
#' breaks its run instead of being bridged by whatever row happens to come next.
#'
#' @param env_dat an `sf` POINT object
#' @param exceeds logical, one per row
#' @param min_steps shortest run kept
#' @param max_gap non-exceeding steps that may be bridged
#' @return a list with `id` (0 outside an event) and `duration`
#' @keywords internal
event_runs <- function(env_dat, exceeds, min_steps = 1L, max_gap = 0L) {
  steps <- time_steps(env_dat)
  step_index <- match_step_index(env_dat, steps)
  cells <- location_key(env_dat)

  id <- integer(length(exceeds))
  duration <- integer(length(exceeds))
  next_id <- 0L

  for (cell in unique(cells)) {
    rows <- which(cells == cell)
    order_in_cell <- rows[order(step_index[rows])]
    present <- step_index[order_in_cell]

    flag <- exceeds[order_in_cell]
    # Bridge short interruptions, but only where the steps either side really
    # are adjacent in time.
    if (max_gap > 0) flag <- bridge_gaps(flag, present, max_gap)

    runs <- rle(flag)
    ends <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L

    for (k in seq_along(runs$values)) {
      if (!isTRUE(runs$values[k])) next
      span <- starts[k]:ends[k]
      # Contiguous in time, not merely adjacent in the ordering.
      contiguous <- all(diff(present[span]) == 1)
      if (!contiguous || length(span) < min_steps) next
      next_id <- next_id + 1L
      id[order_in_cell[span]] <- next_id
      duration[order_in_cell[span]] <- length(span)
    }
  }
  list(id = id, duration = duration)
}

#' Join event runs separated by a short interruption
#'
#' @param flag logical vector in time order
#' @param present the step index of each element
#' @param max_gap longest interruption that may be bridged
#' @return the flag vector with short gaps set `TRUE`
#' @keywords internal
bridge_gaps <- function(flag, present, max_gap) {
  runs <- rle(flag)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L

  for (k in seq_along(runs$values)) {
    interior <- k > 1 && k < length(runs$values)
    if (!interior || isTRUE(runs$values[k])) next
    if (runs$lengths[k] > max_gap) next
    span <- starts[k]:ends[k]
    # Only a gap in time, not a gap in the record.
    if (all(diff(present[(starts[k] - 1L):(ends[k] + 1L)]) == 1)) {
      flag[span] <- TRUE
    }
  }
  flag
}

#' Sum a value over each labelled event
#'
#' @param id event labels, 0 outside an event
#' @param values a numeric vector
#' @return a numeric vector of the event total, broadcast to its rows
#' @keywords internal
event_sum <- function(id, values) {
  totals <- tapply(values, id, function(z) sum(z, na.rm = TRUE))
  out <- as.numeric(totals[match(id, names(totals))])
  out[id == 0] <- NA_real_
  out
}

#' Warn when the threshold rests on too few values
#'
#' @param groups the climatological grouping
#' @param percentile the requested percentile
#' @return invisible `NULL`
#' @keywords internal
warn_thin_threshold <- function(groups, percentile) {
  per_group <- table(groups)
  if (length(per_group) == 0) return(invisible(NULL))
  # Below this many values the percentile is just the extreme of the sample.
  needed <- max(3, ceiling(1 / min(percentile, 1 - percentile)))

  if (max(per_group) < needed) {
    warning(
      "The ", format(100 * percentile), "th percentile is being taken over at ",
      "most ", max(per_group), " value(s) per cell-month, so it is close to ",
      "the extreme of the sample\n  rather than a threshold, and a large ",
      "fraction of steps will be flagged.",
      "\n  Hobday et al. use a 30-year baseline. Fetch a longer series before ",
      "reading anything\n  into the events.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
