#' Find eddies as objects, not just as a field
#'
#' [flow_deformation()] returns the Okubo-Weiss parameter as a number per cell.
#' This takes the next step: it groups the connected cells where rotation beats
#' strain into individual eddies, and describes each one. A cell then carries
#' not "how eddy-like is the flow here" but "you are inside an eddy, it turns
#' this way, and it is this big".
#'
#' That distinction matters ecologically, because the two polarities do opposite
#' things. A cyclonic eddy upwells at its core, lifting nutrients and often
#' concentrating plankton; an anticyclonic one downwells, and its core is
#' typically poorer. A covariate that only says "eddy" averages the two together
#' and can easily find nothing.
#'
#' @section How eddies are identified:
#' The Okubo-Weiss criterion of Isern-Fontanet et al. (2003): a cell belongs to
#' an eddy where \eqn{W < -\alpha\sigma_W}, with \eqn{\sigma_W} the standard
#' deviation of \eqn{W} over that time step and \eqn{\alpha} the `threshold`
#' argument. The published value is 0.2, which is the default. It is a
#' *relative* threshold, recomputed per step, so a quiet month still yields
#' eddies — the criterion asks which parts of this flow are most rotational, not
#' whether the flow is energetic in absolute terms.
#'
#' Connected regions are then labelled, and those smaller than `min_cells` are
#' discarded as noise. On a 1/12-degree product an eddy of oceanographic
#' interest spans many cells; a two-cell patch is more likely to be a numerical
#' artefact of the differencing than a feature.
#'
#' @section What each measure is:
#' * **in_eddy** is 1 inside an identified eddy and 0 outside.
#' * **polarity** is +1 for cyclonic rotation and -1 for anticyclonic, from the
#'   sign of the eddy's mean vorticity. In the northern hemisphere cyclonic is
#'   counter-clockwise and upwelling at the core.
#' * **radius** is the equivalent radius of the eddy, \eqn{\sqrt{A/\pi}} for its
#'   area \eqn{A}, in `per` units. Real eddies are not discs, so this is a size,
#'   not a shape.
#'
#' Cells outside any eddy get 0 for `in_eddy` and `NA` for the others, because
#' they have no eddy to describe.
#'
#' @section What this is not:
#' Detection per time step, with no identity between steps. Nothing here tracks
#' an eddy through its life, so there is no age, no lifespan and no propagation
#' speed, and the same physical eddy carries unrelated labels in consecutive
#' steps. `eddy_id` is therefore deliberately not returned: it would invite
#' exactly that mistake.
#'
#' The Okubo-Weiss criterion is also known to be permissive in strongly strained
#' flow and sensitive to the smoothness of the velocity field. It finds
#' rotation-dominated regions, which is a good working definition of an eddy and
#' not the same thing as a closed-contour eddy from altimetry.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessCopernicus()`
#' @param u name of the eastward velocity column, in m/s
#' @param v name of the northward velocity column, in m/s
#' @param threshold multiple of the standard deviation of Okubo-Weiss below
#'   which a cell counts as rotational. Positive; 0.2 follows the literature
#' @param min_cells smallest connected region kept, in cells
#' @param measures any of `"in_eddy"`, `"polarity"` and `"radius"`
#' @param per distance unit for the radius, `"km"` or `"m"`
#' @param suffix appended to each measure to name its column
#' @return `env_dat` with one column per requested measure
#' @references
#' Isern-Fontanet J, Garcia-Ladona E, Font J (2003). Identification of marine
#' eddies from altimetric maps. *Journal of Atmospheric and Oceanic Technology*
#' **20**(5), 772-778.
#' \doi{10.1175/1520-0426(2003)20<772:IOMEFA>2.0.CO;2}
#' @examples
#' \dontrun{
#' env <- detect_eddies(env)
#'
#' # Cyclonic cores upwell, so they are the ones worth separating out.
#' upwelling <- env$in_eddy == 1 & env$polarity > 0
#' }
#' @seealso [flow_deformation()], [distance_to_eddy()], [eke()]
#' @export
detect_eddies <- function(env_dat, u = "UO", v = "VO", threshold = 0.2,
                          min_cells = 4L,
                          measures = c("in_eddy", "polarity", "radius"),
                          per = c("km", "m"), suffix = "") {
  per <- match.arg(per)
  measures <- match.arg(measures, several.ok = TRUE)
  resolve_vars(env_dat, c(u, v), kind = "spatial")
  check_eddy_args(threshold, min_cells)

  per_time_step(env_dat, c(u, v), function(rast) {
    found <- eddy_patches(rast, u, v, threshold, min_cells, per)

    layers <- list(
      in_eddy = terra::ifel(is.na(found$patches), 0, 1),
      polarity = found$polarity,
      radius = found$radius
    )[measures]

    out <- do.call(c, unname(layers))
    names(out) <- paste0(measures, suffix)
    out
  })
}

#' Distance to the nearest eddy
#'
#' The sibling of [distance_to_front()] and [distance_to_isobath()]. Being
#' outside an eddy but close to one is a different place from being far from
#' any, and a binary inside/outside flag cannot express it: a station 5 km from
#' a rotating core and one 300 km away both score zero.
#'
#' `polarity` narrows it to one kind. Because cyclonic and anticyclonic eddies
#' do opposite things to the water column, distance to the nearest *cyclonic*
#' eddy is often the covariate that carries a signal where distance to the
#' nearest eddy of any sort does not.
#'
#' @inheritParams detect_eddies
#' @param polarity `"any"`, `"cyclonic"`, or `"anticyclonic"`
#' @param name name for the new column
#' @return `env_dat` with a distance column, in `per` units. A cell inside a
#'   qualifying eddy has a distance of 0. A step containing no qualifying eddy
#'   returns `NA` rather than a distance to nothing
#' @examples
#' \dontrun{
#' env <- distance_to_eddy(env, polarity = "cyclonic")
#' }
#' @seealso [detect_eddies()], [distance_to_front()]
#' @export
distance_to_eddy <- function(env_dat, u = "UO", v = "VO", threshold = 0.2,
                             min_cells = 4L,
                             polarity = c("any", "cyclonic", "anticyclonic"),
                             per = c("km", "m"), name = NULL) {
  polarity <- match.arg(polarity)
  per <- match.arg(per)
  resolve_vars(env_dat, c(u, v), kind = "spatial")
  check_eddy_args(threshold, min_cells)

  wanted <- switch(polarity, any = 0, cyclonic = 1, anticyclonic = -1)
  scale <- if (per == "km") 1000 else 1

  result <- per_time_step(env_dat, c(u, v), function(rast) {
    found <- eddy_patches(rast, u, v, threshold, min_cells, per)

    cores <- if (wanted == 0) {
      terra::ifel(is.na(found$patches), NA, 1)
    } else {
      terra::ifel(found$polarity == wanted, 1, NA)
    }

    # An empty mask has nothing to measure from, and terra::distance() would
    # return the distance to the edge of the raster rather than say so.
    if (all(is.na(terra::values(cores)))) {
      out <- terra::ifel(is.na(cores), NA, NA)
    } else {
      out <- terra::distance(cores) / scale
    }
    names(out) <- "eddy_distance__tmp"
    out
  })

  target <- name %||% switch(polarity,
                             any = "eddy_dist",
                             cyclonic = "cyclonic_eddy_dist",
                             anticyclonic = "anticyclonic_eddy_dist")
  result[[target]] <- result[["eddy_distance__tmp"]]
  result[["eddy_distance__tmp"]] <- NULL
  result
}

#' Label rotation-dominated patches within one time step
#'
#' @param rast a `SpatRaster` carrying the velocity components
#' @param u,v names of the velocity layers
#' @param threshold multiple of sd(W) below which a cell is rotational
#' @param min_cells smallest patch kept
#' @param per `"km"` or `"m"`, for the radius
#' @return a list of `SpatRaster`s: `patches`, `polarity`, `radius`
#' @keywords internal
eddy_patches <- function(rast, u, v, threshold, min_cells, per = "km") {
  du <- gradient_layers(rast[[u]], per = "m")
  dv <- gradient_layers(rast[[v]], per = "m")
  vorticity <- dv[[2]] - du[[3]]
  normal <- du[[2]] - dv[[3]]
  shear <- dv[[2]] + du[[3]]
  okubo <- normal^2 + shear^2 - vorticity^2

  blank <- terra::ifel(is.na(okubo), NA, NA)
  empty <- list(patches = blank, polarity = blank, radius = blank)

  spread <- terra::global(okubo, "sd", na.rm = TRUE)[1, 1]
  if (!is.finite(spread) || spread == 0) return(empty)

  rotational <- terra::ifel(okubo < -threshold * spread, 1, NA)
  if (all(is.na(terra::values(rotational)))) return(empty)

  patches <- terra::patches(rotational, directions = 8, allowGaps = FALSE)

  counts <- terra::freq(patches)
  if (!any(counts$count >= min_cells)) return(empty)
  # `%in%` on a SpatRaster gives a plain logical vector rather than a raster,
  # so small patches are removed by substituting their ids away instead.
  small <- counts$value[counts$count < min_cells]
  if (length(small) > 0) {
    patches <- terra::subst(patches, from = small, to = NA)
  }

  # Mean vorticity gives the sense of rotation; total area gives the size.
  spin <- terra::zonal(vorticity, patches, "mean", na.rm = TRUE)
  area <- terra::zonal(terra::cellSize(patches, unit = per), patches, "sum",
                       na.rm = TRUE)
  ids <- spin[[1]]

  list(
    patches = patches,
    polarity = terra::subst(patches, from = ids, to = sign(spin[[2]])),
    radius = terra::subst(patches, from = ids, to = sqrt(area[[2]] / pi))
  )
}

#' Validate the arguments shared by the eddy functions
#'
#' @param threshold multiple of sd(W)
#' @param min_cells smallest patch kept
#' @return invisible `NULL`
#' @keywords internal
check_eddy_args <- function(threshold, min_cells) {
  if (!is.numeric(threshold) || length(threshold) != 1 ||
      !is.finite(threshold) || threshold <= 0) {
    stop("`threshold` must be a single positive number, a multiple of the ",
         "standard deviation of\n  Okubo-Weiss. The sign is applied for you: ",
         "0.2 means W < -0.2 sigma.", call. = FALSE)
  }
  if (!is.numeric(min_cells) || length(min_cells) != 1 || is.na(min_cells) ||
      min_cells < 1 || min_cells != round(min_cells)) {
    stop("`min_cells` must be a single whole number of at least 1.",
         call. = FALSE)
  }
  invisible(NULL)
}
