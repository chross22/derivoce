#' Volume transport across a section
#'
#' Integrates the component of the flow normal to a line, giving transport per
#' unit depth in m^2/s. Multiply by a layer thickness for a volume flux.
#'
#' The units match what the moored-array literature reports, but the quantity
#' does not: those estimates integrate over the full depth of a section, and this
#' integrates one model level along its length. Expect magnitudes an order of
#' magnitude or more apart, and read the output as relative variability rather
#' than as a flux to be compared with a published figure.
#'
#' Unlike most of this package the result is **one number per time step**,
#' broadcast to every row. It describes the section, not the cell, in the same
#' way a climate index describes a basin. A horizontal gradient of it is
#' therefore identically zero, and [horizontal_gradient()] will say so.
#'
#' @section Which way is positive:
#' The normal points to the **right of the direction of travel** from `from` to
#' `to`. Walking the section from `from` to `to`, flow crossing left to right
#' counts positive. Swap the endpoints to reverse the sign.
#'
#' Check this once against a field whose direction you know. A sign error here is
#' invisible: the magnitudes stay plausible and only the interpretation inverts.
#'
#' @section Sampling and gaps:
#' The section is divided into equal segments and the flow is sampled at each
#' midpoint, bilinearly. The default spacing is half a grid cell, so the section
#' is not under-resolved relative to the data it is drawn on.
#'
#' Sample points on land or outside the domain have no velocity. They are
#' dropped rather than counted as zero flow, since zero would understate the
#' transport while looking like a measurement. If fewer than `min_coverage` of
#' the points survive, the step returns `NA`: a transport integrated over half a
#' section is not that section's transport.
#'
#' @section What this is not:
#' A surface-velocity field integrated along a line is a **proxy for**, not a
#' measurement of, the depth-integrated transport that a mooring array gives.
#' Baroclinic structure means the surface can flow one way while the deep channel
#' flows the other, which is exactly the situation in the Northeast Channel. Read
#' the output as an index of variability rather than as a flux.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return, on a
#'   regular lon/lat grid, with eastward and northward velocity columns
#' @param from,to endpoints of the section, each `c(longitude, latitude)`
#' @param u name of the eastward velocity column, in m/s
#' @param v name of the northward velocity column, in m/s
#' @param spacing sample spacing along the section, in km. `NULL` uses half a
#'   grid cell.
#' @param min_coverage fraction of sample points that must carry a velocity for
#'   the step to return a value
#' @param name name for the new column
#' @return `env_dat` with a transport column added, in m^2/s
#' @references
#' Ramp SR, Schlitz RJ, Wright WR (1985). The deep flow through the Northeast
#' Channel, Gulf of Maine. *Journal of Physical Oceanography* **15**(12),
#' 1790-1808.
#' @examples
#' \dontrun{
#' env <- datamatch::accessCopernicus(vars = c("UO", "VO"), ...)
#'
#' # An arbitrary section
#' env <- section_transport(env, from = c(-66.5, 43.3), to = c(-65.6, 42.6))
#'
#' # The named indices, whose geometry is fixed
#' env <- scotian_shelf_inflow(env)
#' env <- northeast_channel_inflow(env)
#' }
#' @seealso [scotian_shelf_inflow()], [northeast_channel_inflow()],
#'   [derived_indices()]
#' @export
section_transport <- function(env_dat, from, to, u = "UO", v = "VO",
                              spacing = NULL, min_coverage = 0.5,
                              name = "transport") {
  resolve_vars(env_dat, c(u, v))
  check_endpoint(from, "from")
  check_endpoint(to, "to")
  if (isTRUE(all.equal(as.numeric(from), as.numeric(to)))) {
    stop("`from` and `to` are the same point, so the section has no length.",
         call. = FALSE)
  }
  if (min_coverage < 0 || min_coverage > 1) {
    stop("min_coverage must be between 0 and 1.", call. = FALSE)
  }

  geometry <- section_geometry(env_dat, from, to, spacing)
  steps <- time_steps(env_dat)
  result <- rep(NA_real_, nrow(env_dat))
  covered <- 0

  for (i in seq_len(nrow(steps))) {
    rows <- step_rows(env_dat, steps[i, ])
    velocity <- rasterize_step(env_dat[rows, ], c(u, v))
    sampled <- as.matrix(terra::extract(velocity, geometry$points,
                                        method = "bilinear"))

    # Normal component at each sample point, in m/s.
    normal <- sampled[, 1] * geometry$normal[1] + sampled[, 2] * geometry$normal[2]
    usable <- !is.na(normal)

    if (mean(usable) >= min_coverage && any(usable)) {
      result[rows] <- sum(normal[usable]) * geometry$ds
      covered <- covered + 1
    }
  }

  if (covered == 0) {
    warning("section_transport() returned no values at all: every point is NA.",
            "\n  No time step had velocities at ", format(min_coverage * 100),
            "% of the section's sample points. The usual cause is a section ",
            "that lies outside the fetched bounding box, or mostly over land.",
            "\n  Check the endpoints against the data's extent. Coordinates are ",
            "c(longitude, latitude), so a transposed pair lands somewhere ",
            "unintended rather than erroring.", call. = FALSE)
  }

  env_dat[[name]] <- result
  env_dat
}

#' Scotian Shelf inflow at Cape Sable
#'
#' Transport of the Nova Scotia Current across a line off Cape Sable, where cold
#' fresh Scotian Shelf Water rounds the cape and enters the Gulf of Maine.
#' Positive values are **into the Gulf**. Informally this is the "Scotian Shelf
#' crossover": the partition of shelf water between entering the Gulf and
#' continuing past it.
#'
#' Scotian Shelf Water supplies more than half of the Gulf of Maine's freshwater
#' budget, and its inflow drives Gulf salinity at seasonal and interannual scales
#' (Feng et al. 2016; Wang et al. 2022). The inflow is also the fresh, cold,
#' nutrient-poor counterpart to the slope water entering through the Northeast
#' Channel, and the two alternate episodically (Townsend et al. 2015), which is
#' why [northeast_channel_inflow()] is a separate index rather than a variant of
#' this one.
#'
#' @section Why the section is fixed:
#' An index named for a place is defined by that place. A "crossover" computed
#' somewhere else is a different quantity that happens to share a function, so
#' the endpoints are not an argument. Use [section_transport()] for any other
#' line, and [scotian_shelf_inflow_section()] to see or plot the geometry used
#' here.
#'
#' The endpoints run from Cape Sable southwest across the shelf toward the head
#' of the Northeast Channel. They are **approximate**, chosen to cut the current
#' where it rounds the cape rather than to reproduce a particular published
#' section, and no published index is being replicated here.
#'
#' @inheritParams section_transport
#' @param name name for the new column
#' @return `env_dat` with a `scotian_inflow` column, in m^2/s, positive into the
#'   Gulf of Maine
#' @references
#' Feng H, Vandemark D, Wilkin J (2016). Gulf of Maine salinity variation and its
#' correlation with upstream Scotian Shelf currents at seasonal and interannual
#' time scales. *Journal of Geophysical Research: Oceans* **121**.
#' \doi{10.1002/2016JC012337}
#'
#' Wang et al. (2022). Freshwater transport in the Scotian Shelf and its impacts
#' on the Gulf of Maine salinity. *Journal of Geophysical Research: Oceans*
#' **127**. \doi{10.1029/2021JC017663}
#' @examples
#' \dontrun{
#' env <- scotian_shelf_inflow(env)
#' }
#' @seealso [northeast_channel_inflow()], [section_transport()]
#' @export
scotian_shelf_inflow <- function(env_dat, u = "UO", v = "VO", spacing = NULL,
                                 min_coverage = 0.5, name = "scotian_inflow") {
  section <- scotian_shelf_inflow_section()
  section_transport(env_dat, from = section$from, to = section$to,
                    u = u, v = v, spacing = spacing,
                    min_coverage = min_coverage, name = name)
}

#' Northeast Channel inflow
#'
#' Transport across the Northeast Channel, between Georges Bank and Browns Bank.
#' Positive values are **into the Gulf of Maine**.
#'
#' A separate index from [scotian_shelf_inflow()], not a variant of it. The
#' Northeast Channel is the deep connection through which slope water enters the
#' Gulf, and it is the Gulf's main source of dissolved inorganic nutrients (Ramp
#' et al. 1985; Townsend et al. 2015). Scotian Shelf inflow at Cape Sable is
#' shallow, fresh, and nutrient-poor. The two vary independently and episodically,
#' so combining them into one number would hide the contrast worth having.
#'
#' @section A caution specific to this section:
#' Ramp et al. (1985) measured the *deep* flow here, below 75 m, and found
#' persistent in-channel transport. Surface velocities over the same section can
#' run the other way, because the channel is strongly baroclinic. An index built
#' from surface currents is therefore not comparable in sign or magnitude to the
#' moored estimates, and is best read as relative variability. If depth-resolved
#' velocities are available, fetch them at channel depth and pass them here
#' instead.
#'
#' @section The inflow regime is not stationary:
#' Slope water entering here is modulated by Gulf Stream warm-core rings. Du et
#' al. (2022) found interannual ring activity off the Gulf tracks bottom salinity
#' in the Northeast Channel, through coastal-trapped waves excited when a ring
#' meets the shelf edge.
#'
#' Silver et al. (2023) then showed the forcing itself changed: ring formation
#' nearly doubled after 2000, from about 18 a year to 33, and salinity-maximum
#' intrusions onto the Northeast Shelf quadrupled, with 72% of observed
#' intrusions coinciding with a ring offshore.
#'
#' So a long record of this index spans two regimes rather than one, and a model
#' fitted across the break may be averaging over a change in the mechanism. This
#' matters most for anything interannual: check whether a relationship holds
#' before and after 2000 separately, rather than assuming it is stable.
#'
#' It also bears on interpretation. A high value here can mean the same
#' circulation carrying more slope water because more rings are present, rather
#' than the circulation itself having strengthened.
#'
#' @inheritParams section_transport
#' @param name name for the new column
#' @return `env_dat` with a `channel_inflow` column, in m^2/s, positive into the
#'   Gulf of Maine
#' @references
#' Ramp SR, Schlitz RJ, Wright WR (1985). The deep flow through the Northeast
#' Channel, Gulf of Maine. *Journal of Physical Oceanography* **15**(12),
#' 1790-1808.
#'
#' Townsend DW, Pettigrew NR, Thomas MA, Neary MG, McGillicuddy DJ, O'Donnell J
#' (2015). Water masses and nutrient sources to the Gulf of Maine. *Journal of
#' Marine Research* **73**, 93-122.
#'
#' Du J, Zhang WG, Li Y (2022). Impact of Gulf Stream warm-core rings on slope
#' water intrusion into the Gulf of Maine. *Journal of Physical Oceanography*
#' **52**(8). \doi{10.1175/JPO-D-21-0288.1}
#'
#' Silver A, Gangopadhyay A, Gawarkiewicz G, Fratantoni P, Clark J (2023).
#' Increased Gulf Stream warm core ring formations contributes to an observed
#' increase in salinity maximum intrusions on the Northeast Shelf. *Scientific
#' Reports* **13**, 7538. \doi{10.1038/s41598-023-34494-0}
#' @examples
#' \dontrun{
#' env <- northeast_channel_inflow(env)
#' }
#' @seealso [scotian_shelf_inflow()], [section_transport()]
#' @export
northeast_channel_inflow <- function(env_dat, u = "UO", v = "VO",
                                     spacing = NULL, min_coverage = 0.5,
                                     name = "channel_inflow") {
  section <- northeast_channel_section()
  section_transport(env_dat, from = section$from, to = section$to,
                    u = u, v = v, spacing = spacing,
                    min_coverage = min_coverage, name = name)
}

#' Endpoints of the named sections
#'
#' Exposed so the geometry behind each index can be inspected, plotted, or
#' adjusted, rather than sitting as a constant inside a function body. Neither
#' reproduces a specific published section.
#'
#' @section How these were placed:
#' Not by eye. Both were measured against GLORYS monthly surface velocities for
#' January to April 2010 on two diagnostics:
#'
#' \itemize{
#'   \item **capture fraction**, `|mean flow . n|` over `|mean flow|`, which is 1
#'     when the current crosses the section squarely and 0 when it runs along it
#'   \item **endpoint ratio**, the larger endpoint normal velocity over the peak
#'     along the section, which is near zero when the section spans the current
#'     rather than cutting through it
#' }
#'
#' An earlier pair chosen from a map scored 0.65 and 0.27 on capture. The
#' Northeast Channel one ran nearly along the channel axis, so its transport was
#' the difference between two opposing halves. The current pair score 0.99 and
#' 0.93, with Channel endpoints at about 4% of the peak.
#'
#' The Northeast Channel line was then cross-checked against ETOPO depth, since
#' bathymetry is what defines that channel: along it, depth runs roughly 120 m,
#' 250 m, 80 m, so it starts on one bank, crosses the deep water, and ends on the
#' other.
#'
#' Both normals point westerly, into the Gulf, which was verified explicitly
#' rather than inferred from the endpoint order.
#'
#' `docs/methods.md` records the full derivation, and
#' `docs/section-placement-diagnostics.R` re-runs it on any `uo`/`vo` extract.
#'
#' @return a list with `from` and `to`, each `c(longitude, latitude)`
#' @examples
#' scotian_shelf_inflow_section()
#' northeast_channel_section()
#' @seealso [section_transport()] for an arbitrary line
#' @export
scotian_shelf_inflow_section <- function() {
  # A near-meridional line at about 66 W, from the Nova Scotia shore south
  # across the passage between Cape Sable and Browns Bank. Travelling from ->
  # to, the right-hand normal points west, which is the way water rounding the
  # cape enters the Gulf, so inflow is positive.
  #
  # Kept after a 60-month re-check, in preference to an orientation that scored
  # higher on capture. That one aligned better with the total flow but nearly
  # cancelled the winter inflow, dropping the seasonal transport from about
  # +1,400 to +160 m^2/s and leaving only 45% of winters positive. Capturing the
  # signal matters more than capturing the flow.
  list(from = c(-66.15, 43.54), to = c(-65.95, 42.76))
}

#' @rdname scotian_shelf_inflow_section
#' @export
northeast_channel_section <- function() {
  # Across the channel at about 66.4 W, from the Browns Bank side down to
  # Georges Bank. Depth along it runs roughly 75 m, 250 m, 75 m, so it starts on
  # a bank, crosses the deep channel, and ends on the other bank, which is what
  # makes it a crossing rather than a chord.
  #
  # Placed on 60 months rather than one season. An earlier version tuned on four
  # months of 2010 scored 0.93 on capture for that window and only 0.59 across
  # five winters; this scores 0.84 there and 0.80 over all 60 months, with a
  # winter transport nearly four times larger, so it is measuring the inflow
  # rather than a residual of it.
  list(from = c(-66.19, 42.72), to = c(-66.61, 41.88))
}

#' Validate a section endpoint
#'
#' @param point a candidate `c(longitude, latitude)`
#' @param label argument name, for the message
#' @return invisibly `TRUE`
#' @keywords internal
check_endpoint <- function(point, label) {
  if (!is.numeric(point) || length(point) != 2 || any(!is.finite(point))) {
    stop("`", label, "` must be a numeric c(longitude, latitude).", call. = FALSE)
  }
  if (abs(point[2]) > 90) {
    stop("`", label, "` has a latitude outside [-90, 90]. The order is ",
         "c(longitude, latitude), the opposite of the usual spoken order.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Sample points, unit normal, and segment length for a section
#'
#' Geometry is done on an equirectangular plane about the section's mean
#' latitude, the same metric the Lagrangian code uses, so a degree of longitude
#' and a degree of latitude are comparable lengths before any normal is taken.
#'
#' @param env_dat an `sf` POINT object, used for the grid spacing
#' @param from,to endpoints, `c(longitude, latitude)`
#' @param spacing sample spacing in km, or `NULL` for half a grid cell
#' @return list with `points` (matrix of lon/lat), `normal` (unit vector), and
#'   `ds` (segment length in metres)
#' @keywords internal
section_geometry <- function(env_dat, from, to, spacing = NULL) {
  metres_per_degree <- 111320
  reference_lat <- mean(c(from[2], to[2]))
  cos_ref <- cos(reference_lat * pi / 180)

  dx <- (to[1] - from[1]) * metres_per_degree * cos_ref
  dy <- (to[2] - from[2]) * metres_per_degree
  length_m <- sqrt(dx^2 + dy^2)

  # Right-hand normal to the direction of travel.
  normal <- c(dy, -dx) / length_m

  step_m <- (spacing %||% (grid_spacing_km(env_dat) / 2)) * 1000
  n <- max(2, ceiling(length_m / step_m))

  # Midpoint rule: sample at segment centres, so no sample sits exactly on an
  # endpoint, where the section may already be off the grid.
  fraction <- (seq_len(n) - 0.5) / n
  points <- cbind(from[1] + fraction * (to[1] - from[1]),
                  from[2] + fraction * (to[2] - from[2]))

  list(points = points, normal = normal, ds = length_m / n)
}

#' Smaller grid cell dimension, in km
#'
#' @param env_dat an `sf` POINT object on a regular grid
#' @return cell size in km
#' @keywords internal
grid_spacing_km <- function(env_dat) {
  xy <- sf::st_coordinates(env_dat)
  lon <- sort(unique(xy[, 1]))
  lat <- sort(unique(xy[, 2]))
  if (length(lon) < 2 || length(lat) < 2) {
    stop("The grid has fewer than two distinct longitudes or latitudes.",
         call. = FALSE)
  }

  metres_per_degree <- 111320
  cos_ref <- cos(mean(lat) * pi / 180)
  min((lon[2] - lon[1]) * metres_per_degree * cos_ref,
      (lat[2] - lat[1]) * metres_per_degree) / 1000
}
