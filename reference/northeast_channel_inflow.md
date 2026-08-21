# Northeast Channel inflow

Transport across the Northeast Channel, between Georges Bank and Browns
Bank. Positive values are **into the Gulf of Maine**.

## Usage

``` r
northeast_channel_inflow(
  env_dat,
  u = "UO",
  v = "VO",
  spacing = NULL,
  min_coverage = 0.5,
  name = "channel_inflow"
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

`env_dat` with a `channel_inflow` column, in m^2/s, positive into the
Gulf of Maine

## Details

A separate index from
[`scotian_shelf_inflow()`](https://camilleross.org/derivoce/reference/scotian_shelf_inflow.md),
not a variant of it. The Northeast Channel is the deep connection
through which slope water enters the Gulf, and it is the Gulf's main
source of dissolved inorganic nutrients (Ramp et al. 1985; Townsend et
al. 2015). Scotian Shelf inflow at Cape Sable is shallow, fresh, and
nutrient-poor. The two vary independently and episodically, so combining
them into one number would hide the contrast worth having.

## A caution specific to this section

Ramp et al. (1985) measured the *deep* flow here, below 75 m, and found
persistent in-channel transport. Surface velocities over the same
section can run the other way, because the channel is strongly
baroclinic. An index built from surface currents is therefore not
comparable in sign or magnitude to the moored estimates, and is best
read as relative variability. If depth-resolved velocities are
available, fetch them at channel depth and pass them here instead.

## The inflow regime is not stationary

Slope water entering here is modulated by Gulf Stream warm-core rings.
Du et al. (2022) found interannual ring activity off the Gulf tracks
bottom salinity in the Northeast Channel, through coastal-trapped waves
excited when a ring meets the shelf edge.

Silver et al. (2023) then showed the forcing itself changed: ring
formation nearly doubled after 2000, from about 18 a year to 33, and
salinity-maximum intrusions onto the Northeast Shelf quadrupled, with
72% of observed intrusions coinciding with a ring offshore.

So a long record of this index spans two regimes rather than one, and a
model fitted across the break may be averaging over a change in the
mechanism. This matters most for anything interannual: check whether a
relationship holds before and after 2000 separately, rather than
assuming it is stable.

It also bears on interpretation. A high value here can mean the same
circulation carrying more slope water because more rings are present,
rather than the circulation itself having strengthened.

## References

Ramp SR, Schlitz RJ, Wright WR (1985). The deep flow through the
Northeast Channel, Gulf of Maine. *Journal of Physical Oceanography*
**15**(12), 1790-1808.

Townsend DW, Pettigrew NR, Thomas MA, Neary MG, McGillicuddy DJ,
O'Donnell J (2015). Water masses and nutrient sources to the Gulf of
Maine. *Journal of Marine Research* **73**, 93-122.

Du J, Zhang WG, Li Y (2022). Impact of Gulf Stream warm-core rings on
slope water intrusion into the Gulf of Maine. *Journal of Physical
Oceanography* **52**(8).
[doi:10.1175/JPO-D-21-0288.1](https://doi.org/10.1175/JPO-D-21-0288.1)

Silver A, Gangopadhyay A, Gawarkiewicz G, Fratantoni P, Clark J (2023).
Increased Gulf Stream warm core ring formations contributes to an
observed increase in salinity maximum intrusions on the Northeast Shelf.
*Scientific Reports* **13**, 7538.
[doi:10.1038/s41598-023-34494-0](https://doi.org/10.1038/s41598-023-34494-0)

## See also

[`scotian_shelf_inflow()`](https://camilleross.org/derivoce/reference/scotian_shelf_inflow.md),
[`section_transport()`](https://camilleross.org/derivoce/reference/section_transport.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- northeast_channel_inflow(env)
} # }
```
