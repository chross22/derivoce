# Marine heatwaves and cold spells

Flags periods when a cell was unusually warm — or unusually cold — for
the time of year, and measures how unusual and for how long. Follows the
definition of Hobday et al. (2016): a threshold set by a high percentile
of the cell's own seasonally varying climatology, an event being a run
of consecutive steps above it, and an intensity measured from the
climatology rather than from the threshold.

## Usage

``` r
marine_heatwave(
  env_dat,
  var = "SST",
  percentile = 0.9,
  min_steps = 1L,
  max_gap = 0L,
  direction = c("warm", "cold"),
  measures = c("event", "intensity", "category", "duration", "cumulative"),
  prefix = NULL
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- var:

  the covariate to examine, typically sea surface temperature

- percentile:

  threshold percentile, between 0 and 1. Hobday et al. use 0.9 for
  heatwaves; the cold-spell equivalent is 0.1, which `direction` applies
  for you

- min_steps:

  shortest run of consecutive steps counted as an event. Use 5 on daily
  data for the published definition

- max_gap:

  runs separated by at most this many non-exceeding steps are joined
  into one event. Use 2 on daily data for the published definition

- direction:

  `"warm"` for heatwaves, `"cold"` for cold spells

- measures:

  which of `"event"`, `"intensity"`, `"category"`, `"duration"` and
  `"cumulative"` to add

- prefix:

  stem for the new column names. Defaults to `"mhw"`, or `"mcs"` when
  `direction = "cold"`

## Value

`env_dat` with one column per requested measure

## Details

The Gulf of Maine is one of the fastest-warming shelf seas on record, so
whether an animal was there in an ordinary summer or a heatwave summer
is often a more useful covariate than the temperature itself.

## What each measure is

- **event** is `TRUE` where the cell is in an event that met
  `min_steps`.

- **intensity** is the departure from the climatological mean, in the
  units of `var`, and is `NA` outside an event. Measured from the
  climatology, not from the threshold, so it is comparable between cells
  whose thresholds differ.

- **category** is 1 to 4 — moderate, strong, severe, extreme — from how
  many multiples of the threshold's excess over the climatology the
  value reaches.

- **duration** is the length of the event in time steps, repeated on
  every row belonging to it.

- **cumulative** is the summed intensity over the whole event, which
  distinguishes a long mild event from a short violent one.

## Time steps, not days

Hobday et al. define an event as five or more consecutive **days**. This
works in whatever step the data is in, because datamatch's access
functions serve monthly, daily and sub-daily archives alike. On daily
data `min_steps = 5` and `max_gap = 2` reproduce the published
definition. On monthly data they do not — five consecutive months is a
far rarer and larger thing — so the defaults here are deliberately
permissive and the choice is left to you.

Runs are counted over consecutive entries in the sequence of time steps
**present in the data**, which has two consequences worth knowing.

A cell absent from a step that other cells have breaks that cell's
event: it cannot be shown to have stayed warm through a step it has no
value for.

A step missing from the record *entirely* is a different matter. It is
not in the sequence at all, so the steps either side of it are adjacent
and an event runs straight through. On a record with holes in it that
can join what were two events, and the duration will be in steps you
have rather than in elapsed time. Check the series is complete before
reading duration as a period.

## Why it needs a long series

A 90th percentile taken over three values is the largest of the three,
and every third step will then be a heatwave. The published definition
uses a 30-year baseline. This warns when the groups behind the threshold
are too thin for the percentile to mean anything, but it cannot tell you
that a 10-year baseline is short for your purpose.

Because the climatology is computed from the data given, a warming trend
within the series raises the threshold and later heatwaves are measured
against a warmer baseline. That is a choice, not an oversight: it makes
events relative to recent conditions. If you want a fixed baseline,
compute the threshold on a subset and apply it.

## References

Hobday AJ, Alexander LV, Perkins SE, Smale DA, Straub SC, Oliver ECJ,
Benthuysen JA, Burrows MT, Donat MG, Feng M, Holbrook NJ, Moore PJ,
Scannell HA, Sen Gupta A, Wernberg T (2016). A hierarchical approach to
defining marine heatwaves. *Progress in Oceanography* **141**, 227-238.
[doi:10.1016/j.pocean.2015.12.014](https://doi.org/10.1016/j.pocean.2015.12.014)

Hobday AJ, Oliver ECJ, Sen Gupta A, Benthuysen JA, Burrows MT, Donat MG,
Holbrook NJ, Moore PJ, Thomsen MS, Wernberg T, Smale DA (2018).
Categorizing and naming marine heatwaves. *Oceanography* **31**(2),
162-173.
[doi:10.5670/oceanog.2018.205](https://doi.org/10.5670/oceanog.2018.205)

## See also

[`cell_anomaly()`](https://camilleross.org/derivoce/reference/cell_anomaly.md),
[`box_anomaly()`](https://camilleross.org/derivoce/reference/box_anomaly.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- marine_heatwave(env, "SST")

# The published definition, on daily data.
env <- marine_heatwave(env, "SST", min_steps = 5, max_gap = 2)

# Cold spells instead.
env <- marine_heatwave(env, "SST", direction = "cold")
} # }
```
