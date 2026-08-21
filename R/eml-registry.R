#' What each fixed-name derived column means
#'
#' Columns this package creates under a name of its own, rather than by
#' suffixing a source column. Kept beside the suffix registry so that adding a
#' derivation and forgetting to describe it is one obvious omission rather than
#' a silent one.
#'
#' @return a named list of `definition`, `unit` and `number`
#' @keywords internal
eml_fixed_registry <- function() {
  list(
    speed = list(
      definition = "Current speed, the magnitude of the horizontal velocity vector, sqrt(u^2 + v^2).",
      unit = "metersPerSecond"),
    EKE = list(
      definition = "Eddy kinetic energy, 0.5 * (u'^2 + v'^2), where the anomalies are departures from a mean whose definition is set by the `reference` argument.",
      unit = "metersSquaredPerSecondSquared"),
    shore_dist = list(
      definition = "Great-circle distance to the nearest coastline, from Natural Earth.",
      unit = "kilometer"),
    transport = list(
      definition = "Volume transport across a section, positive in the direction of the section normal.",
      unit = "cubicMetersPerSecond"),
    scotian_inflow = list(
      definition = "Volume transport of Scotian Shelf Water across the Cape Sable section into the Gulf of Maine.",
      unit = "cubicMetersPerSecond"),
    channel_inflow = list(
      definition = "Volume transport across the Northeast Channel section.",
      unit = "cubicMetersPerSecond"),
    egom_salinity = list(
      definition = "Eastern Gulf of Maine surface salinity anomaly, an index of winter Scotian Shelf inflow after Grodsky et al. (2025).",
      unit = "practicalSalinityUnit"),
    sigma_theta = list(
      definition = "Potential density anomaly (sigma-theta): density at one atmosphere from the UNESCO (1983) equation of state, minus 1000.",
      unit = "kilogramPerCubicMeter"),
    density = list(
      definition = "Potential density at one atmosphere, from the UNESCO (1983) equation of state.",
      unit = "kilogramPerCubicMeter"),
    N2 = list(
      definition = "Buoyancy frequency squared, -(g/rho0) * d(rho)/dz, between two depths. Positive where the column is stably stratified.",
      unit = "perSecondSquared"),
    N = list(
      definition = "Buoyancy frequency between two depths.",
      unit = "perSecond"),
    eady_growth = list(
      definition = "Eady growth rate of baroclinic instability, 0.31 * |f| * |dU/dz| / N, in the maximum-growth form of Lindzen and Farrell (1980) for the instability of Eady (1949).",
      unit = "perDay"),
    vorticity = list(
      definition = "Relative vorticity, dv/dx - du/dy. Positive is counter-clockwise.",
      unit = "perSecond"),
    divergence = list(
      definition = "Horizontal divergence, du/dx + dv/dy. Negative is convergence.",
      unit = "perSecond"),
    normal_strain = list(
      definition = "Normal component of strain, du/dx - dv/dy.",
      unit = "perSecond"),
    shear_strain = list(
      definition = "Shear component of strain, dv/dx + du/dy.",
      unit = "perSecond"),
    strain_rate = list(
      definition = "Total strain rate, the magnitude of the two strain components.",
      unit = "perSecond"),
    okubo_weiss = list(
      definition = "Okubo-Weiss parameter, Sn^2 + Ss^2 - zeta^2. Negative where rotation dominates, which is the interior of a coherent eddy.",
      unit = "perSecondSquared"),
    rossby = list(
      definition = "Rossby number, relative vorticity divided by the Coriolis parameter.",
      unit = "dimensionless"),
    in_eddy = list(
      definition = "Whether the cell lies inside an identified eddy: 1 inside, 0 outside.",
      unit = "dimensionless", number = "integer"),
    polarity = list(
      definition = "Sense of rotation of the eddy containing this cell: +1 cyclonic, -1 anticyclonic.",
      unit = "dimensionless", number = "integer"),
    radius = list(
      definition = "Equivalent radius of the eddy containing this cell, sqrt(area/pi).",
      unit = "kilometer"),
    eddy_dist = list(
      definition = "Distance to the nearest identified eddy, zero inside one.",
      unit = "kilometer"),
    cyclonic_eddy_dist = list(
      definition = "Distance to the nearest cyclonic eddy, zero inside one.",
      unit = "kilometer"),
    anticyclonic_eddy_dist = list(
      definition = "Distance to the nearest anticyclonic eddy, zero inside one.",
      unit = "kilometer"),
    forward_ftle = list(
      definition = "Forward finite-time Lyapunov exponent, marking repelling structures and transport barriers.",
      unit = "perDay"),
    backward_ftle = list(
      definition = "Backward finite-time Lyapunov exponent, marking attracting structures where material accumulates.",
      unit = "perDay"),
    forward_fsle = list(
      definition = "Forward finite-size Lyapunov exponent.",
      unit = "perDay"),
    backward_fsle = list(
      definition = "Backward finite-size Lyapunov exponent.",
      unit = "perDay"),
    forward_residence = list(
      definition = "Time water starting in this cell remains inside the region before leaving. Right-censored at the tracking limit.",
      unit = "nominalDay"),
    backward_residence = list(
      definition = "Time water arriving in this cell had already spent inside the region. Right-censored at the tracking limit.",
      unit = "nominalDay"),
    mhw_event = list(
      definition = "Whether the cell is within a marine heatwave or cold spell, after Hobday et al. (2016): 1 within, 0 outside.",
      unit = "dimensionless", number = "integer"),
    mhw_category = list(
      definition = "Marine heatwave category after Hobday et al. (2018): 1 moderate, 2 strong, 3 severe, 4 extreme.",
      unit = "dimensionless", number = "integer"),
    mhw_duration = list(
      definition = "Length of the marine heatwave containing this cell, in time steps of the input series.",
      unit = "number", number = "whole")
  )
}

#' What each suffixed derived column means
#'
#' Matched against the end of a column name, longest patterns first so that
#' `_front_dist` is not mistaken for a rolling statistic. Each entry supplies a
#' function giving the unit, because most derived units are built from the
#' source column's own unit.
#'
#' @return a list of pattern entries
#' @keywords internal
eml_suffix_registry <- function() {
  list(
    list(pattern = "_front_dist$",
         definition = "Distance to the nearest front, found by thresholding the horizontal gradient after Belkin and O'Reilly (2009).",
         unit = function(source, units) "kilometer"),
    list(pattern = "_front_freq$",
         definition = "Fraction of time steps in which this cell was itself frontal.",
         unit = function(source, units) "dimensionless"),
    list(pattern = "_box_anom$",
         definition = "Departure of the box mean from its reference, one value per time step broadcast to every row.",
         unit = function(source, units) derived_unit(source, units)),
    list(pattern = "_grad_x$",
         definition = "Eastward component of the horizontal gradient, by central differences.",
         unit = function(source, units) derived_unit(source, units, "PerKilometer")),
    list(pattern = "_grad_y$",
         definition = "Northward component of the horizontal gradient, by central differences.",
         unit = function(source, units) derived_unit(source, units, "PerKilometer")),
    list(pattern = "_grad$",
         definition = "Magnitude of the horizontal gradient, by central differences on the lon/lat lattice, per unit distance rather than per degree.",
         unit = function(source, units) derived_unit(source, units, "PerKilometer")),
    list(pattern = "_vgrad$",
         definition = "Vertical difference between two levels, a stratification index.",
         unit = function(source, units) NA_character_),
    list(pattern = "_tgrad$",
         definition = "Rate of change between consecutive time steps.",
         unit = function(source, units) derived_unit(source, units)),
    list(pattern = "_int$",
         definition = "Covariate accumulated over preceding time steps.",
         unit = function(source, units) derived_unit(source, units)),
    list(pattern = "_seasonal$",
         definition = "Repeating seasonal component, estimated jointly with the trend and centred on zero.",
         unit = function(source, units) derived_unit(source, units)),
    list(pattern = "_residual$",
         definition = "What remains once the mean, the trend and the seasonal cycle are removed.",
         unit = function(source, units) derived_unit(source, units)),
    list(pattern = "_trend$",
         definition = "Fitted long-term component, centred on zero.",
         unit = function(source, units) derived_unit(source, units)),
    list(pattern = "_slope$",
         definition = "Rate of long-term change, one value per cell.",
         unit = function(source, units) derived_unit(source, units, "PerYear")),
    list(pattern = "_anom$",
         definition = "Departure from this cell's own mean, over the whole record or within calendar month.",
         unit = function(source, units) derived_unit(source, units)),
    list(pattern = "_z$",
         definition = "Departure from this cell's own mean, divided by that cell's standard deviation.",
         unit = function(source, units) "dimensionless"),
    list(pattern = "_frac$",
         definition = "Fraction of the water present originating from the named endmember, by temperature-salinity mixing.",
         unit = function(source, units) "dimensionless"),
    list(pattern = "_lag[0-9]+(step|day|month|year)?$",
         definition = "Covariate lagged by a fixed number of positions or of calendar units.",
         unit = function(source, units) derived_unit(source, units)),
    list(pattern = "_(mean|sd|min|max|sum|median|range)[0-9]+(day|month|year)?$",
         definition = "Summary of a trailing window ending at, and including, the current time step.",
         unit = function(source, units) derived_unit(source, units))
  )
}

#' Declarations for units outside the EML dictionary
#'
#' EML validates against a fixed dictionary of 195 units. These are the ones
#' this package produces that are not in it, in the shape a `unitList` needs.
#'
#' @return a data frame of custom unit declarations
#' @keywords internal
eml_unit_registry <- function() {
  data.frame(
    id = c("perSecond", "perSecondSquared", "perDay",
           "metersSquaredPerSecondSquared", "practicalSalinityUnit",
           "celsiusPerKilometer", "celsiusPerYear", "sverdrup",
           "millimolesPerCubicMeter", "milligramsPerSquareMeterPerDay",
           "milligramsPerCubicMeterPerDay", "wattsPerSquareMeter"),
    unitType = c("frequency", "unknown", "frequency", "unknown", "dimensionless",
                 "unknown", "unknown", "volumetricRate", "amountOfSubstance",
                 "unknown", "unknown", "powerPerArea"),
    parentSI = c("hertz", NA, "hertz", NA, NA, NA, NA,
                 "cubicMetersPerSecond", "molePerCubicMeter", NA, NA, "watt"),
    multiplierToSI = c(1, NA, 1 / 86400, NA, NA, NA, NA,
                       1e6, 0.001, NA, NA, 1),
    description = c(
      "Reciprocal seconds, as for vorticity, divergence, strain rate and buoyancy frequency.",
      "Reciprocal seconds squared, as for the Okubo-Weiss parameter and the square of the buoyancy frequency.",
      "Reciprocal days, as for Lyapunov exponents and the Eady growth rate.",
      "Metres squared per second squared, the unit of kinetic energy per unit mass.",
      "Practical Salinity Unit, dimensionless on the Practical Salinity Scale 1978.",
      "Degrees Celsius per kilometre, as for a horizontal temperature gradient.",
      "Degrees Celsius per year, as for a fitted temperature trend.",
      "Sverdrup, 10^6 cubic metres per second, the conventional unit of ocean volume transport.",
      "Millimoles per cubic metre, as datamatch serves nitrate, phosphate and dissolved oxygen.",
      "Milligrams per square metre per day, as datamatch serves satellite primary production.",
      "Milligrams per cubic metre per day, as datamatch serves modelled net primary production.",
      "Watts per square metre, as FVCOM serves shortwave radiation and net surface heat flux."
    ),
    stringsAsFactors = FALSE
  )
}

#' EML units for the variables datamatch serves
#'
#' A translation of datamatch's own unit strings into EML unit ids, so that
#' [eml_attributes()] can resolve the source columns of a normal workflow
#' without being told what each one holds. Derived units are built on top of
#' these, so a temperature gradient becomes degrees per kilometre without the
#' caller naming either part.
#'
#' Hardcoded rather than read from datamatch at runtime, because datamatch is
#' deliberately not a dependency of this package. The cost is that it can fall
#' behind, so a test checks it against the live catalogue whenever datamatch is
#' installed, and reports any variable this table has not been told about.
#'
#' @return a named character vector, variable name to EML unit id
#' @keywords internal
datamatch_units <- function() {
  c(
    # Copernicus physics
    SST = "celsius", BOTT = "celsius",
    SSS = "practicalSalinityUnit", BOTS = "practicalSalinityUnit",
    UO = "metersPerSecond", VO = "metersPerSecond",
    SSH = "meter", MLD = "meter", SIC = "dimensionless",
    # Surface forcing. Wind stress is N/m^2, which is a pascal exactly.
    WSPD = "metersPerSecond", UWND = "metersPerSecond",
    VWND = "metersPerSecond",
    TAUX = "pascal", TAUY = "pascal", TAU = "pascal",
    # Copernicus biogeochemistry
    CHL = "milligramsPerCubicMeter", CHL_MODEL = "milligramsPerCubicMeter",
    DIATO = "milligramsPerCubicMeter", DINO = "milligramsPerCubicMeter",
    PP = "milligramsPerSquareMeterPerDay",
    NPP_MODEL = "milligramsPerCubicMeterPerDay",
    NO3 = "millimolesPerCubicMeter", PO4 = "millimolesPerCubicMeter",
    O2 = "millimolesPerCubicMeter", PH = "dimensionless",
    # Seafloor terrain, from attach_bathymetry()
    DEPTH = "meter", SLOPE = "degree", ASPECT = "degree", TPI = "meter",
    # Climate indices, from attach_climate_index(). The oscillation indices are
    # standardized anomalies; AMOC is a transport in Sverdrups.
    NAO = "dimensionless", AO = "dimensionless", AMO = "dimensionless",
    PDO = "dimensionless", LCR = "dimensionless", AMOC = "sverdrup",
    # Satellite wind observation counts, from accessCCMP()
    NOBS = "number",
    # HYCOM adds bottom velocities beside the bottom temperature and salinity
    UO_BOTTOM = "metersPerSecond", VO_BOTTOM = "metersPerSecond",
    # FVCOM adds depth-averaged velocities and surface fluxes
    UBAR = "metersPerSecond", VBAR = "metersPerSecond",
    SWRAD = "wattsPerSquareMeter", NHF = "wattsPerSquareMeter"
  )
}
