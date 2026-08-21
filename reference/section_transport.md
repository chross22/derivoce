# Volume transport across a section

Integrates the component of the flow normal to a line, giving transport
per unit depth in m^2/s. Multiply by a layer thickness for a volume
flux.

## Usage

``` r
section_transport(
  env_dat,
  from,
  to,
  u = "UO",
  v = "VO",
  spacing = NULL,
  min_coverage = 0.5,
  name = "transport"
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return, on a regular lon/lat grid, with
  eastward and northward velocity columns

- from, to:

  endpoints of the section, each `c(longitude, latitude)`

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

`env_dat` with a transport column added, in m^2/s

## Details

The units match what the moored-array literature reports, but the
quantity does not: those estimates integrate over the full depth of a
section, and this integrates one model level along its length. Expect
magnitudes an order of magnitude or more apart, and read the output as
relative variability rather than as a flux to be compared with a
published figure.

Unlike most of this package the result is **one number per time step**,
broadcast to every row. It describes the section, not the cell, in the
same way a climate index describes a basin. A horizontal gradient of it
is therefore identically zero, and
[`horizontal_gradient()`](https://camilleross.org/derivoce/reference/horizontal_gradient.md)
will say so.

## Which way is positive

The normal points to the **right of the direction of travel** from
`from` to `to`. Walking the section from `from` to `to`, flow crossing
left to right counts positive. Swap the endpoints to reverse the sign.

Check this once against a field whose direction you know. A sign error
here is invisible: the magnitudes stay plausible and only the
interpretation inverts.

## Sampling and gaps

The section is divided into equal segments and the flow is sampled at
each midpoint, bilinearly. The default spacing is half a grid cell, so
the section is not under-resolved relative to the data it is drawn on.

Sample points on land or outside the domain have no velocity. They are
dropped rather than counted as zero flow, since zero would understate
the transport while looking like a measurement. If fewer than
`min_coverage` of the points survive, the step returns `NA`: a transport
integrated over half a section is not that section's transport.

## What this is not

A surface-velocity field integrated along a line is a **proxy for**, not
a measurement of, the depth-integrated transport that a mooring array
gives. Baroclinic structure means the surface can flow one way while the
deep channel flows the other, which is exactly the situation in the
Northeast Channel. Read the output as an index of variability rather
than as a flux.

## References

Ramp SR, Schlitz RJ, Wright WR (1985). The deep flow through the
Northeast Channel, Gulf of Maine. *Journal of Physical Oceanography*
**15**(12), 1790-1808.

## See also

[`scotian_shelf_inflow()`](https://camilleross.org/derivoce/reference/scotian_shelf_inflow.md),
[`northeast_channel_inflow()`](https://camilleross.org/derivoce/reference/northeast_channel_inflow.md),
[`derived_indices()`](https://camilleross.org/derivoce/reference/derived_indices.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- datamatch::accessCopernicus(vars = c("UO", "VO"), ...)

# An arbitrary section
env <- section_transport(env, from = c(-66.5, 43.3), to = c(-65.6, 42.6))

# The named indices, whose geometry is fixed
env <- scotian_shelf_inflow(env)
env <- northeast_channel_inflow(env)
} # }
```
