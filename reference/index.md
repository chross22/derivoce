# Package index

## Gradients and fronts

A front is where a field changes sharply, not where it is large. These
compute the change and then locate it.

- [`distance_to_contour()`](https://camilleross.org/derivoce/reference/distance_to_contour.md)
  : Distance to a contour of a covariate
- [`distance_to_front()`](https://camilleross.org/derivoce/reference/distance_to_front.md)
  : Distance to the nearest front
- [`distance_to_isobath()`](https://camilleross.org/derivoce/reference/distance_to_isobath.md)
  : Distance to isobaths
- [`distance_to_shore()`](https://camilleross.org/derivoce/reference/distance_to_shore.md)
  : Distance to the nearest shoreline
- [`front_frequency()`](https://camilleross.org/derivoce/reference/front_frequency.md)
  : How often a place is frontal
- [`horizontal_gradient()`](https://camilleross.org/derivoce/reference/horizontal_gradient.md)
  : Horizontal spatial gradient of a covariate
- [`temporal_gradient()`](https://camilleross.org/derivoce/reference/temporal_gradient.md)
  : Temporal gradient of a covariate
- [`vertical_gradient()`](https://camilleross.org/derivoce/reference/vertical_gradient.md)
  : Vertical temperature gradient

## Flow and mixing

Quantities derived from velocity fields: how fast, how deformed, how
strongly stratified, and how long water stays.

- [`buoyancy_frequency()`](https://camilleross.org/derivoce/reference/buoyancy_frequency.md)
  : Buoyancy frequency between two depths
- [`current_speed()`](https://camilleross.org/derivoce/reference/current_speed.md)
  : Current speed from velocity components
- [`eady_growth_rate()`](https://camilleross.org/derivoce/reference/eady_growth_rate.md)
  : Eady growth rate of baroclinic instability
- [`eke()`](https://camilleross.org/derivoce/reference/eke.md) : Eddy
  kinetic energy
- [`flow_deformation()`](https://camilleross.org/derivoce/reference/flow_deformation.md)
  : Velocity gradient diagnostics: vorticity, divergence, strain,
  Okubo-Weiss
- [`potential_density()`](https://camilleross.org/derivoce/reference/potential_density.md)
  : Potential density of seawater
- [`residence_time()`](https://camilleross.org/derivoce/reference/residence_time.md)
  : How long water stays in a place
- [`water_mass_fraction()`](https://camilleross.org/derivoce/reference/water_mass_fraction.md)
  : Water-mass fraction from two endmembers

## Coherent structures

Eddies, and the Lyapunov exponents that reveal the transport barriers
between them.

- [`detect_eddies()`](https://camilleross.org/derivoce/reference/detect_eddies.md)
  : Find eddies as objects, not just as a field
- [`distance_to_eddy()`](https://camilleross.org/derivoce/reference/distance_to_eddy.md)
  : Distance to the nearest eddy
- [`fsle()`](https://camilleross.org/derivoce/reference/fsle.md) :
  Finite-Size Lyapunov Exponent
- [`ftle()`](https://camilleross.org/derivoce/reference/ftle.md) :
  Finite-Time Lyapunov Exponent

## Time treatments

The same field asked about over a different window: lagged by calendar
units, integrated, rolled, or expressed as an anomaly.

- [`box_anomaly()`](https://camilleross.org/derivoce/reference/box_anomaly.md)
  : Anomaly of a covariate averaged over a box
- [`cell_anomaly()`](https://camilleross.org/derivoce/reference/cell_anomaly.md)
  : Anomaly of a covariate against each cell's own history
- [`decompose_covariate()`](https://camilleross.org/derivoce/reference/decompose_covariate.md)
  : Split a covariate into trend, seasonal cycle, and what is left
- [`integrate_covariate()`](https://camilleross.org/derivoce/reference/integrate_covariate.md)
  : Time-integrated covariate
- [`lag_covariate()`](https://camilleross.org/derivoce/reference/lag_covariate.md)
  : Lagged covariate values
- [`marine_heatwave()`](https://camilleross.org/derivoce/reference/marine_heatwave.md)
  : Marine heatwaves and cold spells
- [`rolling_covariate()`](https://camilleross.org/derivoce/reference/rolling_covariate.md)
  : Rolling summary of a covariate over a trailing window

## Named regions and sections

Standing definitions for the Gulf of Maine features this work returns
to, so a section is computed the same way every time.

- [`eastern_gom_box()`](https://camilleross.org/derivoce/reference/eastern_gom_box.md)
  : The eastern Gulf of Maine box
- [`eastern_gom_salinity()`](https://camilleross.org/derivoce/reference/eastern_gom_salinity.md)
  : Eastern Gulf of Maine salinity index
- [`northeast_channel_inflow()`](https://camilleross.org/derivoce/reference/northeast_channel_inflow.md)
  : Northeast Channel inflow
- [`scotian_shelf_inflow()`](https://camilleross.org/derivoce/reference/scotian_shelf_inflow.md)
  : Scotian Shelf inflow at Cape Sable
- [`scotian_shelf_inflow_section()`](https://camilleross.org/derivoce/reference/scotian_shelf_inflow_section.md)
  [`northeast_channel_section()`](https://camilleross.org/derivoce/reference/scotian_shelf_inflow_section.md)
  : Endpoints of the named sections
- [`section_transport()`](https://camilleross.org/derivoce/reference/section_transport.md)
  : Volume transport across a section

## Indices and metadata

Series assembled from the above, and the EML pieces that travel with
them.

- [`derived_indices()`](https://camilleross.org/derivoce/reference/derived_indices.md)
  [`print(`*`<derivoce_index_table>`*`)`](https://camilleross.org/derivoce/reference/derived_indices.md)
  : Regional indices this package can derive
- [`eml_attributes()`](https://camilleross.org/derivoce/reference/eml_attributes.md)
  : Describe derived columns as EML metadata
- [`eml_col_classes()`](https://camilleross.org/derivoce/reference/eml_col_classes.md)
  : Column classes for EML, in the order the attributes are given
- [`eml_custom_units()`](https://camilleross.org/derivoce/reference/eml_custom_units.md)
  : Declarations for the units EML does not know
- [`index_series()`](https://camilleross.org/derivoce/reference/index_series.md)
  : Pull a per-time-step index back out as a series
