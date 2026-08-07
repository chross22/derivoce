suppressPackageStartupMessages({library(sf); library(terra)})
for (f in list.files("/Users/camille/GitHub/derivoce/R", full.names = TRUE)) source(f)

nc <- ncdf4::nc_open("/tmp/gom_wide.nc")
lon <- nc$dim$longitude$vals; lat <- nc$dim$latitude$vals
tv <- nc$dim$time$vals
origin <- as.POSIXct(sub("^hours since ", "", ncdf4::ncatt_get(nc,"time","units")$value), tz="UTC")
dates <- as.Date(origin + tv*3600)
uo <- ncdf4::ncvar_get(nc,"uo"); vo <- ncdf4::ncvar_get(nc,"vo")
ncdf4::nc_close(nc)

env <- st_as_sf(do.call(rbind, lapply(seq_along(dates), function(i) {
  g <- expand.grid(x=lon,y=lat); g$UO <- as.numeric(uo[,,i]); g$VO <- as.numeric(vo[,,i])
  g$YEAR <- as.integer(format(dates[i],"%Y")); g$MONTH <- as.integer(format(dates[i],"%m")); g$DAY <- 1L; g
})), coords=c("x","y"), crs=4326)

steps <- time_steps(env)
rasters <- lapply(seq_len(nrow(steps)), function(i)
  rasterize_step(env[step_rows(env, steps[i,]), ], c("UO","VO")))

score <- function(from, to) {
  geo <- try(section_geometry(env, from, to, spacing = 2), silent = TRUE)
  if (inherits(geo,"try-error")) return(NULL)
  caps <- c(); ends <- c(); nets <- c()
  for (r in rasters) {
    s <- as.matrix(terra::extract(r, geo$points, method="bilinear"))
    if (mean(!is.na(s[,1])) < 0.95) return(NULL)          # must stay over water
    n <- s[,1]*geo$normal[1] + s[,2]*geo$normal[2]
    mu <- c(mean(s[,1],na.rm=TRUE), mean(s[,2],na.rm=TRUE))
    caps <- c(caps, abs(sum(mu*geo$normal))/sqrt(sum(mu^2)))
    ends <- c(ends, max(abs(n[1]), abs(n[length(n)])) / max(abs(n), na.rm=TRUE))
    nets <- c(nets, sum(n, na.rm=TRUE)*geo$ds)
  }
  if (!all(is.finite(caps)) || !all(is.finite(ends)) || !all(is.finite(nets))) return(NULL)
  data.frame(cap=mean(caps), endr=mean(ends), net=mean(nets),
             sign_ok = all(nets > 0), len=max(dim(geo$points)[1])*geo$ds/1000)
}

# Candidate sections: centre, orientation (deg from east, of the SECTION line),
# half-length in degrees of latitude.
search <- function(cx, cy, angles, halves, label) {
  best <- NULL
  for (a in angles) for (h in halves) {
    cosr <- cos(cy*pi/180)
    dy <- h*sin(a*pi/180); dx <- h*cos(a*pi/180)/cosr
    from <- c(cx-dx, cy-dy); to <- c(cx+dx, cy+dy)
    s <- score(from, to)
    if (is.null(s)) next
    # Want: high capture, low endpoint flow, correct sign
    s$obj <- s$cap - 0.5*s$endr
    if (!s$sign_ok) { from2 <- to; to <- from; from <- from2; s2 <- score(from,to)
                      if (is.null(s2) || !s2$sign_ok) next; s <- s2; s$obj <- s$cap - 0.5*s$endr }
    s$from1 <- from[1]; s$from2 <- from[2]; s$to1 <- to[1]; s$to2 <- to[2]; s$angle <- a
    if (!is.finite(s$obj)) next
    if (is.null(best) || s$obj > best$obj) best <- s
  }
  cat(sprintf("\n%s  best placement\n", label))
  cat(sprintf("  from c(%.2f, %.2f)  to c(%.2f, %.2f)   length %.0f km, angle %.0f deg\n",
      best$from1, best$from2, best$to1, best$to2, best$len, best$angle))
  cat(sprintf("  capture %.2f   endpoint/peak %.2f   mean net %+.0f m^2/s\n",
      best$cap, best$endr, best$net))
  best
}

cs <- search(-66.05, 43.15, seq(0, 170, by = 10), seq(0.30, 0.75, by = 0.05), "CAPE SABLE")
nec <- search(-66.40, 42.30, seq(0, 170, by = 10), seq(0.30, 0.60, by = 0.05), "NORTHEAST CHANNEL")
saveRDS(list(cs=cs, nec=nec), "/tmp/best.rds")
