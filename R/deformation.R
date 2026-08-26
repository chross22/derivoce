#' Velocity gradient diagnostics: vorticity, divergence, strain, Okubo-Weiss
#'
#' Derives the local structure of a flow from the four horizontal derivatives of
#' its velocity components. Together these say whether a parcel of water is
#' spinning, converging, being pulled apart, or some combination.
#'
#' All are instantaneous and local: they describe the flow in one cell at one
#' time step, and need no history. That places them between [eke()], which needs
#' a series to form an anomaly, and [ftle()] and [fsle()], which advect particles
#' through many steps. Where those two answer "how variable is this place" and
#' "where do trajectories separate", these answer "what is the flow doing here,
#' right now".
#'
#' @section The measures:
#' Writing the derivatives of eastward `u` and northward `v` with respect to
#' eastward `x` and northward `y`:
#'
#' * **vorticity**, \eqn{\zeta = \partial v/\partial x - \partial u/\partial y},
#'   is rotation. Positive is counter-clockwise (cyclonic in the northern
#'   hemisphere), which in an eddy means upwelling at its core.
#' * **divergence**, \eqn{\partial u/\partial x + \partial v/\partial y}, is
#'   spreading. Negative is convergence, where surface water piles up and sinks,
#'   and where buoyant particles and floating material accumulate.
#' * **normal_strain**, \eqn{\partial u/\partial x - \partial v/\partial y}, and
#'   **shear_strain**, \eqn{\partial v/\partial x + \partial u/\partial y},
#'   stretch a parcel along the axes and diagonally.
#' * **strain_rate** is their magnitude, the total rate of deformation.
#' * **okubo_weiss**, \eqn{W = S_n^2 + S_s^2 - \zeta^2}, compares the two.
#'   Negative means rotation wins, which is the interior of a coherent eddy.
#'   Positive means strain wins, which is the filaments between eddies where
#'   water is drawn out into long thin structures.
#' * **rossby** is \eqn{\zeta/f}, vorticity as a fraction of planetary rotation.
#'   It is the dimensionless version, so it can be compared between latitudes,
#'   and values approaching 1 mean the flow is fast enough for the usual
#'   geostrophic balance to be breaking down.
#'
#' @section What it cannot tell you:
#' These are single-time-step quantities, so a large value means the flow is
#' deforming now, not that it has been or will be. A coherent eddy that persists
#' for weeks and a momentary filament can carry the same Okubo-Weiss value.
#' Persistence is what [ftle()] and [fsle()] measure, at much greater cost.
#'
#' The derivatives are central differences, so like every gradient in this
#' package the outermost ring of cells is `NA`, and the values are only as
#' smooth as the grid. On a coarse product a real eddy smaller than a few cells
#' will not appear at all.
#'
#' @param env_dat an `sf` POINT object with one row per location and time step,
#'   as datamatch's access functions return
#' @param u name of the eastward velocity column, in m/s
#' @param v name of the northward velocity column, in m/s
#' @param measures which diagnostics to add. Any of `"vorticity"`,
#'   `"divergence"`, `"normal_strain"`, `"shear_strain"`, `"strain_rate"`,
#'   `"okubo_weiss"` and `"rossby"`
#' @param suffix appended to each measure to name its column
#' @return `env_dat` with one column per requested measure. Vorticity,
#'   divergence and the strains are in s^-1, Okubo-Weiss in s^-2, and the Rossby
#'   number is dimensionless
#' @references
#' Okubo A (1970). Horizontal dispersion of floatable particles in the vicinity
#' of velocity singularities such as convergences. *Deep-Sea Research and
#' Oceanographic Abstracts* **17**(3), 445-454.
#' \doi{10.1016/0011-7471(70)90059-8}
#'
#' Weiss J (1991). The dynamics of enstrophy transfer in two-dimensional
#' hydrodynamics. *Physica D: Nonlinear Phenomena* **48**(2-3), 273-294.
#' \doi{10.1016/0167-2789(91)90088-Q}
#'
#' Isern-Fontanet J, Garcia-Ladona E, Font J (2003). Identification of marine
#' eddies from altimetric maps. *Journal of Atmospheric and Oceanic Technology*
#' **20**(5), 772-778.
#' \doi{10.1175/1520-0426(2003)20<772:IOMEFA>2.0.CO;2}
#' @examples
#' \dontrun{
#' env <- flow_deformation(env)
#'
#' # Eddy interiors are where rotation beats strain.
#' eddy_core <- env$okubo_weiss < 0
#' }
#' @seealso [eke()], [ftle()], [fsle()], [horizontal_gradient()]
#' @export
flow_deformation <- function(env_dat, u = "UO", v = "VO",
                             measures = c("vorticity", "divergence",
                                          "strain_rate", "okubo_weiss"),
                             suffix = "") {
  known <- c("vorticity", "divergence", "normal_strain", "shear_strain",
             "strain_rate", "okubo_weiss", "rossby")
  measures <- match.arg(measures, choices = known, several.ok = TRUE)
  resolve_vars(env_dat, c(u, v), kind = "spatial")

  per_time_step(env_dat, c(u, v), function(rast) {
    # Metres, not kilometres, so the derivatives come out in s^-1 rather than
    # in a mixed unit that would then have to be carried through Okubo-Weiss.
    du <- gradient_layers(rast[[u]], per = "m")
    dv <- gradient_layers(rast[[v]], per = "m")
    du_dx <- du[[2]]; du_dy <- du[[3]]
    dv_dx <- dv[[2]]; dv_dy <- dv[[3]]

    vorticity <- dv_dx - du_dy
    normal_strain <- du_dx - dv_dy
    shear_strain <- dv_dx + du_dy

    # switch() rather than a list subset, so only the requested measures are
    # evaluated. Building all seven and then dropping six would call coriolis()
    # on every call, and it stops on a projected grid -- which would fail a
    # request for vorticity alone, citing a measure the caller never asked for.
    layers <- lapply(measures, function(measure) {
      switch(measure,
             vorticity = vorticity,
             divergence = du_dx + dv_dy,
             normal_strain = normal_strain,
             shear_strain = shear_strain,
             strain_rate = sqrt(normal_strain^2 + shear_strain^2),
             okubo_weiss = normal_strain^2 + shear_strain^2 - vorticity^2,
             rossby = vorticity / coriolis(rast))
    })

    out <- do.call(c, unname(layers))
    names(out) <- paste0(measures, suffix)
    out
  })
}

#' Coriolis parameter over a raster
#'
#' \eqn{f = 2\Omega\sin(\phi)}, which vanishes at the equator. A grid spanning
#' it would divide by something arbitrarily close to zero, so the band where
#' that happens is returned as `NA` rather than as a very large Rossby number
#' that looks like a real signal.
#'
#' @param rast a `SpatRaster`, used for its geometry only
#' @return a `SpatRaster` of `f`, in s^-1
#' @keywords internal
coriolis <- function(rast) {
  if (!terra::is.lonlat(rast)) {
    stop("The Rossby number needs latitude, so it cannot be computed on a ",
         "projected grid. Drop \"rossby\" from `measures`.", call. = FALSE)
  }
  omega <- 7.2921e-5
  latitude <- terra::init(rast[[1]], "y")
  f <- 2 * omega * sin(latitude * pi / 180)

  # Within about two degrees of the equator f is under 5e-6 s^-1 and the ratio
  # stops meaning anything.
  terra::ifel(abs(latitude) < 2, NA, f)
}
