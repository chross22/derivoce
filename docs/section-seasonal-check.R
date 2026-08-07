suppressPackageStartupMessages({library(sf); library(terra)})
for (f in list.files("/Users/camille/GitHub/derivoce/R", full.names=TRUE)) source(f)

nc <- ncdf4::nc_open("/tmp/gom_5yr.nc")
lon <- nc$dim$longitude$vals; lat <- nc$dim$latitude$vals; tv <- nc$dim$time$vals
origin <- as.POSIXct(sub("^hours since ","",ncdf4::ncatt_get(nc,"time","units")$value), tz="UTC")
dates <- as.Date(origin + tv*3600)
uo <- ncdf4::ncvar_get(nc,"uo"); vo <- ncdf4::ncvar_get(nc,"vo"); ncdf4::nc_close(nc)
cat("months available:", length(dates), "  ", format(min(dates),"%Y-%m"), "to", format(max(dates),"%Y-%m"), "\n\n")

build <- function(idx) st_as_sf(do.call(rbind, lapply(idx, function(i) {
  g <- expand.grid(x=lon,y=lat); g$UO <- as.numeric(uo[,,i]); g$VO <- as.numeric(vo[,,i])
  g$YEAR <- as.integer(format(dates[i],"%Y")); g$MONTH <- as.integer(format(dates[i],"%m")); g$DAY <- 1L; g
})), coords=c("x","y"), crs=4326)

mon <- as.integer(format(dates, "%m"))
seasons <- list(winter = which(mon %in% 1:4), summer = which(mon %in% 6:9))

diag_one <- function(env, s) {
  geo <- section_geometry(env, s$from, s$to, spacing = 2)
  steps <- time_steps(env)
  caps <- ends <- nets <- numeric(0)
  for (i in seq_len(nrow(steps))) {
    r <- rasterize_step(env[step_rows(env, steps[i,]), ], c("UO","VO"))
    smp <- as.matrix(terra::extract(r, geo$points, method="bilinear"))
    n <- smp[,1]*geo$normal[1] + smp[,2]*geo$normal[2]
    mu <- c(mean(smp[,1],na.rm=TRUE), mean(smp[,2],na.rm=TRUE))
    caps <- c(caps, abs(sum(mu*geo$normal))/sqrt(sum(mu^2)))
    ends <- c(ends, max(abs(n[1]),abs(n[length(n)]))/max(abs(n),na.rm=TRUE))
    nets <- c(nets, sum(n,na.rm=TRUE)*geo$ds)
  }
  list(cap=mean(caps), cap_sd=sd(caps), endr=mean(ends),
       net=mean(nets), net_sd=sd(nets), pos=mean(nets > 0))
}

secs <- list("Cape Sable"=scotian_shelf_inflow_section(),
             "Northeast Channel"=northeast_channel_section())

cat(sprintf("%-19s %-8s %6s %7s %8s %11s %7s\n",
            "section","season","capture","(sd)","endpt","net m^2/s","% pos"))
cat(strrep("-", 72), "\n")
for (nm in names(secs)) for (sn in names(seasons)) {
  env <- build(seasons[[sn]])
  d <- diag_one(env, secs[[nm]])
  cat(sprintf("%-19s %-8s %6.2f %7.2f %8.2f %11s %6.0f%%\n",
      nm, sn, d$cap, d$cap_sd, d$endr,
      paste0(format(round(d$net), big.mark=","), " +/-", format(round(d$net_sd), big.mark=",")),
      100*d$pos))
}
