#' @keywords internal
"_PACKAGE"

# Null-coalescing operator. Defined here rather than relying on base R's, which
# only exists from R 4.4 onward.
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Time-step columns that are not covariates
#'
#' @keywords internal
time_columns <- function() c("YEAR", "MONTH", "DAY")

#' Covariate column names in an environmental data object
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`
#' @return character vector of covariate column names
#' @export
covariate_columns <- function(env_dat) {
  setdiff(names(env_dat), c(time_columns(), attr(env_dat, "sf_column")))
}

#' Split environmental data into its time steps
#'
#' @param env_dat an `sf` POINT object with `YEAR`/`MONTH`/`DAY` columns
#' @return a data frame of unique time steps, ordered chronologically
#' @keywords internal
time_steps <- function(env_dat) {
  steps <- unique(sf::st_drop_geometry(env_dat)[time_columns()])
  steps <- steps[order(steps$YEAR, steps$MONTH, steps$DAY), , drop = FALSE]
  rownames(steps) <- NULL
  steps
}

#' Row indices belonging to one time step
#'
#' @param env_dat an `sf` POINT object
#' @param step a one-row data frame with `YEAR`/`MONTH`/`DAY`
#' @return integer vector of row indices
#' @keywords internal
step_rows <- function(env_dat, step) {
  which(env_dat$YEAR == step$YEAR & env_dat$MONTH == step$MONTH &
          env_dat$DAY == step$DAY)
}

#' Rasterize one time step's points onto their own grid
#'
#' Spatial derivatives are only defined on a grid, so the points have to be put
#' back onto one. Gridded ocean products come as a regular lon/lat lattice, and
#' `datamatch::accessEnvDat()` flattens that lattice to points without moving
#' them, so the grid can be recovered exactly from the unique coordinates.
#'
#' Scattered points are rejected rather than interpolated: silently gridding
#' irregular data would produce a gradient field that looks plausible and is
#' mostly interpolation artifact.
#'
#' @param points an `sf` POINT object for a single time step
#' @param vars covariate columns to include as layers
#' @return a `terra::SpatRaster` with one layer per variable
#' @keywords internal
rasterize_step <- function(points, vars) {
  coords <- sf::st_coordinates(points)
  lon <- sort(unique(coords[, 1]))
  lat <- sort(unique(coords[, 2]))

  if (length(lon) < 2 || length(lat) < 2) {
    stop("A time step has fewer than two distinct longitudes or latitudes, ",
         "so no spatial gradient can be computed from it.", call. = FALSE)
  }
  if (!is_regular(lon) || !is_regular(lat)) {
    stop("Points are not on a regular lon/lat grid, so spatial derivatives are ",
         "not well defined. Gridded products (Copernicus and similar) are; ",
         "scattered observations are not.", call. = FALSE)
  }

  values <- as.data.frame(sf::st_drop_geometry(points)[vars])
  frame <- cbind(x = coords[, 1], y = coords[, 2], values)
  terra::rast(frame, type = "xyz", crs = sf::st_crs(points)$wkt)
}

#' Check whether coordinates are evenly spaced
#'
#' @param values sorted unique coordinate values
#' @param tolerance relative tolerance on the spacing
#' @return `TRUE` if the spacing is constant within tolerance
#' @keywords internal
is_regular <- function(values, tolerance = 1e-6) {
  spacing <- diff(values)
  max(abs(spacing - spacing[1])) <= tolerance * abs(spacing[1])
}

#' Apply a raster-valued function to every time step
#'
#' Handles the split/rasterize/compute/join cycle that each spatial derivative
#' shares, so those functions only have to say what to do to one raster.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessEnvDat()`
#' @param vars covariate columns the function needs
#' @param fun a function taking a `SpatRaster` and returning a `SpatRaster`
#'   whose layer names become the new columns
#' @return `env_dat` with the returned layers added as columns
#' @keywords internal
per_time_step <- function(env_dat, vars, fun) {
  steps <- time_steps(env_dat)
  new_columns <- NULL

  for (i in seq_len(nrow(steps))) {
    rows <- step_rows(env_dat, steps[i, ])
    points <- env_dat[rows, ]

    derived <- fun(rasterize_step(points, vars))
    sampled <- terra::extract(derived, sf::st_coordinates(points))

    if (is.null(new_columns)) {
      new_columns <- matrix(NA_real_, nrow = nrow(env_dat), ncol = ncol(sampled),
                            dimnames = list(NULL, names(sampled)))
    }
    new_columns[rows, ] <- as.matrix(sampled)
  }

  for (column in colnames(new_columns)) {
    env_dat[[column]] <- new_columns[, column]
  }
  env_dat
}

#' Check that requested covariates exist
#'
#' @param env_dat an `sf` POINT object
#' @param vars requested covariate names, or `NULL` for all of them
#' @return the resolved covariate names
#' @keywords internal
resolve_vars <- function(env_dat, vars) {
  available <- covariate_columns(env_dat)
  if (is.null(vars)) return(available)

  missing <- setdiff(vars, available)
  if (length(missing) > 0) {
    stop("Covariate(s) not present: ", paste(missing, collapse = ", "),
         "\nAvailable: ", paste(available, collapse = ", "), call. = FALSE)
  }
  vars
}
