# Experiment: can the Labrador Current retroflection index be recomputed from
# Copernicus velocities, to extend it past 2014?
#
# The answer is no, at least from monthly fields. See lcr-extension-experiment.md
# for the finding and the diagnosis. This is kept so the negative result is
# reproducible rather than folklore.
#
# Usage, after fetching monthly GLORYS uo/vo over 40-58N, 70-45W:
#   Rscript docs/lcr-extension-experiment.R <particles> <releases> <days> <step>

# Reproduce the Labrador Current retroflection index from Copernicus velocities.
#
# Jutras et al. (2023) seeded virtual particles weekly across a line on the
# Labrador Shelf, tracked each for three years through GLORYS12V1 with
# OceanParcels, and took the difference between the counts reaching the Labrador
# and Scotian shelves. This is the same idea on monthly GLORYS, using the RK4
# advection already in derivoce.
#
# It is a reimplementation, not their series: different integrator, different
# release cadence, monthly rather than daily fields.

suppressPackageStartupMessages({library(sf); library(terra)})
for (f in list.files("/Users/camille/GitHub/derivoce/R", full.names = TRUE)) source(f)

nc <- ncdf4::nc_open("/tmp/lcr_pilot.nc")
lon <- nc$dim$longitude$vals; lat <- nc$dim$latitude$vals; tv <- nc$dim$time$vals
origin <- as.POSIXct(sub("^hours since ", "",
                         ncdf4::ncatt_get(nc, "time", "units")$value), tz = "UTC")
dates <- as.Date(origin + tv * 3600)
uo <- ncdf4::ncvar_get(nc, "uo"); vo <- ncdf4::ncvar_get(nc, "vo")
ncdf4::nc_close(nc)

# One SpatRaster per month, built the way rasterize_step() does so the float32
# coordinates are accepted.
make_raster <- function(i) {
  r <- terra::rast(nrows = length(lat), ncols = length(lon),
                   xmin = min(lon) - diff(lon[1:2])/2, xmax = max(lon) + diff(lon[1:2])/2,
                   ymin = min(lat) - diff(lat[1:2])/2, ymax = max(lat) + diff(lat[1:2])/2,
                   crs = "EPSG:4326", nlyrs = 2)
  # terra rows run north to south; the NetCDF latitude runs south to north, so
  # the latitude index is reversed. No transpose: the array is [lon, lat], and R
  # flattens column-major, which already gives terra's row order (all longitudes
  # for the northernmost latitude first). Transposing puts every value in the
  # wrong cell, and the field still looks plausible, which is how it survived a
  # first run.
  flip <- rev(seq_along(lat))
  r[[1]] <- as.numeric(uo[, flip, i])
  r[[2]] <- as.numeric(vo[, flip, i])
  names(r) <- c("UO", "VO")
  r
}
velocity <- lapply(seq_along(dates), make_raster)
times <- as.numeric(dates)

# --- trajectories, not just endpoints ---------------------------------------
# advect() returns where a parcel ended. Counting a crossing needs the path, so
# the RK4 step is looped here and every position kept.
track <- function(seeds, start_time, days, step_days) {
  n <- ceiling(days / step_days)
  path <- array(NA_real_, dim = c(nrow(seeds), 2, n + 1))
  path[, , 1] <- seeds
  pos <- seeds
  t <- start_time
  for (k in seq_len(n)) {
    pos <- rk4_step(pos, t, step_days, velocity, times)
    path[, , k + 1] <- pos
    t <- t + step_days
  }
  path
}

# --- did a path ever reach a region? ----------------------------------------
reaches <- function(path, box) {
  x <- path[, 1, ]; y <- path[, 2, ]
  inside <- x >= box$xmin & x <= box$xmax & y >= box$ymin & y <= box$ymax
  inside[is.na(inside)] <- FALSE
  apply(inside, 1, any)
}

# Downstream on the Labrador Shelf, and the Scotian Shelf. Both approximate.
labrador <- list(xmin = -56.0, xmax = -52.0, ymin = 50.0, ymax = 52.0)
scotian  <- list(xmin = -66.0, xmax = -59.0, ymin = 42.5, ymax = 45.5)

release <- function(n) {
  f <- seq(0, 1, length.out = n)
  cbind(-56.7 + f * (-52.0 + 56.7), 53.0 + f * (54.3 - 53.0))
}

args <- commandArgs(trailingOnly = TRUE)
n_particles <- as.integer(args[1] %||% 60)
n_releases  <- as.integer(args[2] %||% 6)
track_days  <- as.numeric(args[3] %||% 1095)
step_days   <- as.numeric(args[4] %||% 1)

seeds <- release(n_particles)
# The western end of the line is at the Labrador coast, so some seeds land on
# land cells and have no velocity. Drop them rather than tracking NAs.
on_water <- !is.na(terra::extract(velocity[[1]][[1]], seeds)[, 1])
seeds <- seeds[on_water, , drop = FALSE]
n_particles <- nrow(seeds)
cat(sprintf("seeds %d over water (of %s), releases %d, %.0f days each, step %.1f d\n",
            n_particles, length(on_water), n_releases, track_days, step_days))

result <- data.frame()
for (i in seq_len(n_releases)) {
  t0 <- Sys.time()
  path <- track(seeds, times[i], track_days, step_days)
  lab <- reaches(path, labrador)
  sco <- reaches(path, scotian)
  lost <- is.na(path[, 1, dim(path)[3]])
  result <- rbind(result, data.frame(
    YEAR = as.integer(format(dates[i], "%Y")),
    MONTH = as.integer(format(dates[i], "%m")),
    n = n_particles, labrador = sum(lab), scotian = sum(sco),
    lost = sum(lost),
    secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)))
  cat(sprintf("  %s  labrador %3d  scotian %3d  lost %3d  (%.0fs)\n",
              format(dates[i], "%Y-%m"), sum(lab), sum(sco), sum(lost),
              tail(result$secs, 1)))
}
result$index <- (result$labrador - result$scotian) / result$n
saveRDS(result, "/tmp/lcr_result.rds")
cat("\nindex range:", round(range(result$index), 3), "\n")
