#' Distance to the nearest shoreline
#'
#' Distance from each point to the nearest coast, in kilometres. Static — it does
#' not vary in time — so it is computed once for the unique locations and shared
#' across every time step.
#'
#' Distance from shore is a broad proxy for several things at once: depth,
#' terrestrial nutrient input, tidal mixing, and larval retention all covary with
#' it. That makes it a useful covariate and a poor explanation — a model leaning
#' on it is telling you *where*, not *why*. This is the covariate the older
#' pipeline had as `dist`.
#'
#' @section Choosing a coastline:
#' `resolution` selects the Natural Earth coastline detail:
#'
#' \itemize{
#'   \item `"medium"` (the default, 1:50m) suits shelf-scale work.
#'   \item `"large"` (1:10m) resolves individual islands and inlets, which
#'     matters in a place like the Gulf of Maine where the coast is deeply
#'     indented, but is slower and needs `rnaturalearthhires`.
#'   \item `"small"` (1:110m) is too coarse for coastal work; it smooths bays
#'     away entirely.
#' }
#'
#' Distances are computed against true geodesic distance on the ellipsoid rather
#' than a projected approximation, so they stay correct across a wide domain.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`
#' @param resolution Natural Earth coastline resolution: `"medium"`, `"large"`,
#'   or `"small"`
#' @param margin degrees of padding around the data when cropping the coastline.
#'   Land just outside the study area can still be the nearest shore, so cropping
#'   tightly to the bounding box would overstate distances at the edges.
#' @param name name for the new column
#' @return `env_dat` with a distance-to-shore column, in km
#' @examples
#' \dontrun{
#' env <- distance_to_shore(env)
#' # Finer coastline, for a deeply indented shore:
#' env <- distance_to_shore(env, resolution = "large")
#' }
#' @export
distance_to_shore <- function(env_dat, resolution = c("medium", "large", "small"),
                              margin = 5, name = "shore_dist") {
  resolution <- match.arg(resolution)

  if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
    stop("The 'rnaturalearth' package is required for distance to shore. ",
         "Install it with install.packages('rnaturalearth').", call. = FALSE)
  }

  coastline <- shoreline(env_dat, resolution = resolution, margin = margin)

  # Distance to shore is static and the grid repeats every time step, so
  # computing per unique location rather than per row saves the whole
  # calculation on all but the first month.
  coords <- sf::st_coordinates(env_dat)
  key <- paste(coords[, 1], coords[, 2], sep = ",")
  unique_index <- !duplicated(key)
  unique_points <- sf::st_geometry(env_dat)[unique_index]

  nearest <- sf::st_nearest_feature(unique_points, coastline)
  distances <- sf::st_distance(unique_points, coastline[nearest, ],
                               by_element = TRUE)
  distances <- as.numeric(distances) / 1000

  env_dat[[name]] <- distances[match(key, key[unique_index])]
  env_dat
}

#' Coastline geometry covering the data's extent
#'
#' @param env_dat an `sf` object
#' @param resolution Natural Earth resolution
#' @param margin degrees of padding around the data
#' @return an `sf` geometry of coastlines
#' @keywords internal
shoreline <- function(env_dat, resolution = "medium", margin = 5) {
  if (resolution == "large" && !requireNamespace("rnaturalearthhires", quietly = TRUE)) {
    stop("resolution = \"large\" needs the 'rnaturalearthhires' package:\n",
         "  install.packages('rnaturalearthhires', ",
         "repos = 'https://ropensci.r-universe.dev')", call. = FALSE)
  }

  coastline <- rnaturalearth::ne_coastline(scale = resolution, returnclass = "sf")
  coastline <- sf::st_transform(coastline, sf::st_crs(env_dat))

  bbox <- sf::st_bbox(env_dat)
  bbox["xmin"] <- bbox["xmin"] - margin
  bbox["xmax"] <- bbox["xmax"] + margin
  bbox["ymin"] <- max(bbox["ymin"] - margin, -90)
  bbox["ymax"] <- min(bbox["ymax"] + margin, 90)

  cropped <- suppressWarnings(sf::st_crop(coastline, bbox))
  if (nrow(cropped) == 0) {
    stop("No coastline within ", margin, " degrees of the data. Increase ",
         "`margin`, or check the data's coordinates.", call. = FALSE)
  }
  sf::st_geometry(cropped)
}
