# Scotian Shelf inflow at Cape Sable

Transport of the Nova Scotia Current across a line off Cape Sable, where
cold fresh Scotian Shelf Water rounds the cape and enters the Gulf of
Maine. Positive values are **into the Gulf**. Informally this is the
"Scotian Shelf crossover": the partition of shelf water between entering
the Gulf and continuing past it.

## Usage

``` r
scotian_shelf_inflow(
  env_dat,
  u = "UO",
  v = "VO",
  spacing = NULL,
  min_coverage = 0.5,
  name = "scotian_inflow"
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return, on a regular lon/lat grid, with
  eastward and northward velocity columns

- u:

  name of the eastward velocity column, in m/s

- v:

  name of the northward velocity column, in m/s

- spacing:

  sample spacing along the section, in km. `NULL` uses half a grid cell.

- min_coverage:

  fraction of sample points that must carry a velocity for the step to
  return a value

- name:

  name for the new column

## Value

`env_dat` with a `scotian_inflow` column, in m^2/s, positive into the
Gulf of Maine

## Details

Scotian Shelf Water supplies more than half of the Gulf of Maine's
freshwater budget, and its inflow drives Gulf salinity at seasonal and
interannual scales (Feng et al. 2016; Wang et al. 2022). The inflow is
also the fresh, cold, nutrient-poor counterpart to the slope water
entering through the Northeast Channel, and the two alternate
episodically (Townsend et al. 2015), which is why
[`northeast_channel_inflow()`](https://camilleross.org/derivoce/reference/northeast_channel_inflow.md)
is a separate index rather than a variant of this one.

## Why the section is fixed

An index named for a place is defined by that place. A "crossover"
computed somewhere else is a different quantity that happens to share a
function, so the endpoints are not an argument. Use
[`section_transport()`](https://camilleross.org/derivoce/reference/section_transport.md)
for any other line, and
[`scotian_shelf_inflow_section()`](https://camilleross.org/derivoce/reference/scotian_shelf_inflow_section.md)
to see or plot the geometry used here.

The endpoints run from Cape Sable southwest across the shelf toward the
head of the Northeast Channel. They are **approximate**, chosen to cut
the current where it rounds the cape rather than to reproduce a
particular published section, and no published index is being replicated
here.

## References

Feng H, Vandemark D, Wilkin J (2016). Gulf of Maine salinity variation
and its correlation with upstream Scotian Shelf currents at seasonal and
interannual time scales. *Journal of Geophysical Research: Oceans*
**121**.
[doi:10.1002/2016JC012337](https://doi.org/10.1002/2016JC012337)

Wang et al. (2022). Freshwater transport in the Scotian Shelf and its
impacts on the Gulf of Maine salinity. *Journal of Geophysical Research:
Oceans* **127**.
[doi:10.1029/2021JC017663](https://doi.org/10.1029/2021JC017663)

## See also

[`northeast_channel_inflow()`](https://camilleross.org/derivoce/reference/northeast_channel_inflow.md),
[`section_transport()`](https://camilleross.org/derivoce/reference/section_transport.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- scotian_shelf_inflow(env)
} # }
```
