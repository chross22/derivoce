#' @keywords internal
"_PACKAGE"

#' Time-step columns that are not covariates
#'
#' @keywords internal
time_columns <- function() c("YEAR", "MONTH", "DAY")

#' Covariate column names in an environmental data object
#'
#' Deliberately internal, and a deliberate duplicate of
#' `datamatch::covariate_columns()`, which is identical. Exporting it would mask
#' datamatch's whenever both packages are attached, for no gain: callers who want
#' it have it from datamatch already. It is not imported from there because
#' datamatch is a suggested rather than a hard dependency — everything in this
#' package operates on the *shape* `accessEnvDat()` returns, not on datamatch
#' itself, so an object of that shape from any source works.
#'
#' `<var>_source` is included and `<var>_depth` is not, matching
#' `datamatch::covariate_columns()`. A source tag travels with the variable it
#' describes and survives resampling as a categorical; a depth is the model
#' level a derived bottom value was taken from, and the mean of two of those is
#' not the depth any value came from. `.datamatch_source` is datamatch's
#' internal per-row tag and is never data.
#'
#' @param env_dat an `sf` POINT object from `datamatch::accessCopernicus()`
#' @return character vector of covariate column names
#' @keywords internal
covariate_columns <- function(env_dat) {
  candidates <- setdiff(names(env_dat),
                        c(time_columns(), attr(env_dat, "sf_column")))
  candidates[!is_bookkeeping_column(candidates)]
}

#' Columns that are bookkeeping rather than data
#'
#' Kept identical to `datamatch::is_bookkeeping_column()`. Duplicated rather
#' than imported for the same reason `covariate_columns()` is: this package
#' works on the shape datamatch returns, not on datamatch.
#'
#' @param names <char> column names to inspect
#' @return <logical> one per name
#' @keywords internal
is_bookkeeping_column <- function(names) {
  grepl("_depth$", names) | names == ".datamatch_source"
}

#' Every column a deposited table should describe
#'
#' Wider than [covariate_columns()], and deliberately so. `<var>_depth` must not
#' be derived from or resampled as though it were a measurement, but an archive
#' still has to say what it is - `datamatch::write_eml()` documents it, with
#' metres as its unit, rather than leaving a column in the deposited table with
#' nothing describing it. Only datamatch's internal per-row tag is left out,
#' since it does not survive `matchData()` to reach a deposit at all.
#'
#' @param env_dat an `sf` POINT object
#' @return character vector of column names
#' @keywords internal
documented_columns <- function(env_dat) {
  candidates <- setdiff(names(env_dat),
                        c(time_columns(), attr(env_dat, "sf_column")))
  candidates[candidates != ".datamatch_source"]
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

  # The raster is built from the lattice's dimensions rather than handed the raw
  # coordinates, because terra::rast(type = "xyz") runs its own regularity check
  # and rejects the same float32 quantisation that is_regular() has just decided
  # to tolerate. Constructing the grid explicitly puts the two in agreement:
  # having established the spacing is uniform, the median is what it is.
  nx <- length(lon)
  ny <- length(lat)
  dx <- stats::median(diff(lon))
  dy <- stats::median(diff(lat))

  grid <- terra::rast(
    nrows = ny, ncols = nx,
    xmin = min(lon) - dx / 2, xmax = max(lon) + dx / 2,
    ymin = min(lat) - dy / 2, ymax = max(lat) + dy / 2,
    crs = sf::st_crs(points)$wkt, nlyrs = length(vars)
  )

  # terra numbers rows from the top, so the northernmost latitude is row 1.
  column <- match(coords[, 1], lon)
  row <- ny - match(coords[, 2], lat) + 1
  cell <- (row - 1) * nx + column

  values <- sf::st_drop_geometry(points)[vars]
  for (i in seq_along(vars)) {
    layer <- rep(NA_real_, nx * ny)
    layer[cell] <- as.numeric(values[[vars[i]]])
    grid[[i]] <- layer
  }
  names(grid) <- vars
  grid
}

#' Check whether coordinates are evenly spaced
#'
#' The tolerance has to clear the quantisation in the coordinates themselves.
#' Copernicus stores longitude and latitude as **float32**, whose spacing near
#' 67 degrees resolves to about 8e-6, so a nominally uniform 1/12-degree grid
#' arrives with spacings varying by around 1e-4 of a cell. That is a property of
#' the file format, not of the grid, and a tolerance tight enough to reject it
#' rejects every real Copernicus download.
#'
#' 1e-3 is loose enough for that and still far tighter than any genuine
#' irregularity. Scattered observations have spacings that differ by order one,
#' and a grid missing a row or column has one interval of double width, so both
#' are still caught by a wide margin.
#'
#' The comparison is against the median spacing rather than the first, so a
#' single odd interval at the start cannot set the reference for everything else.
#'
#' @param values sorted unique coordinate values
#' @param tolerance relative tolerance on the spacing
#' @return `TRUE` if the spacing is constant within tolerance
#' @keywords internal
is_regular <- function(values, tolerance = 1e-3) {
  spacing <- diff(values)
  reference <- stats::median(spacing)
  if (!is.finite(reference) || reference == 0) return(FALSE)
  max(abs(spacing - reference)) <= tolerance * abs(reference)
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

#' Check that requested covariates exist and can be operated on
#'
#' Three checks, deliberately of different severities:
#'
#' A **missing** column is an error: nothing can be done.
#'
#' A **non-numeric** column is an error when named explicitly — a factor cannot be
#' differentiated or summed, so the request cannot be honoured — but is skipped
#' silently when `vars = NULL` swept it up. A caller who did not name
#' `CHL_source` did not mean it, and failing the whole call over a column they
#' never asked for would make the `NULL` default unusable on any object that has
#' been through `datamatch::fill_satellite_gaps()`.
#'
#' A **degenerate** column — static where the operation is temporal, spatially
#' uniform where it is spatial — is only a warning, because the computation is
#' well defined and the caller may have meant it. See [warn_degenerate()].
#'
#' @param env_dat an `sf` POINT object
#' @param vars requested covariate names, or `NULL` for all of them
#' @param kind the kind of operation the caller is about to perform, so the
#'   degeneracy relevant to it can be checked: `"temporal"`, `"spatial"`, or
#'   `"any"` (the default) for operations that are neither, such as a column-wise
#'   difference
#' @return the resolved covariate names
#' @keywords internal
resolve_vars <- function(env_dat, vars, kind = c("any", "spatial", "temporal")) {
  kind <- match.arg(kind)
  available <- covariate_columns(env_dat)
  is_numeric_column <- function(v) is.numeric(env_dat[[v]])

  if (is.null(vars)) {
    vars <- available[vapply(available, is_numeric_column, logical(1))]
    if (length(vars) == 0) {
      stop("No numeric covariate columns to work on.",
           "\nColumns present: ", paste(available, collapse = ", "),
           call. = FALSE)
    }
  } else {
    missing <- setdiff(vars, available)
    if (length(missing) > 0) {
      stop("Covariate(s) not present: ", paste(missing, collapse = ", "),
           "\nAvailable: ", paste(available, collapse = ", "), call. = FALSE)
    }
    non_numeric <- vars[!vapply(vars, is_numeric_column, logical(1))]
    if (length(non_numeric) > 0) {
      stop("Covariate(s) not numeric: ", paste(non_numeric, collapse = ", "),
           "\nThese cannot be differentiated, lagged, or summed. ",
           "datamatch::fill_satellite_gaps() adds a `<var>_source` factor of ",
           "this kind, recording where each value came from.", call. = FALSE)
    }
  }

  warn_degenerate(env_dat, vars, kind)
  vars
}
