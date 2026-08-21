# How often a place is frontal

The fraction of time steps in which a cell's own gradient was sharp
enough to count as a front. Where
[`distance_to_front()`](https://camilleross.org/derivoce/reference/distance_to_front.md)
asks how far the nearest front was at one moment, this asks how reliably
a front sits in this place at all.

## Usage

``` r
front_frequency(
  env_dat,
  var,
  threshold = NULL,
  quantile = 0.9,
  scope = c("record", "step"),
  per = c("km", "m"),
  n = NULL,
  by = c("step", "day", "month", "year"),
  name = NULL
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- var:

  the covariate whose gradient defines a front

- threshold:

  gradient magnitude at or above which a cell is frontal. When `NULL`,
  taken from `quantile`

- quantile:

  quantile of the gradient used as the threshold

- scope:

  `"record"` for one cutoff across the series, `"step"` for a cutoff per
  time step

- per:

  distance unit for the gradient, `"km"` or `"m"`

- n:

  length of the trailing window, or `NULL` for the whole record

- by:

  `"step"`, `"day"`, `"month"` or `"year"`, as in
  [`rolling_covariate()`](https://camilleross.org/derivoce/reference/rolling_covariate.md)

- name:

  name for the new column

## Value

`env_dat` with a frequency column between 0 and 1

## Details

The distinction matters because fronts move. A cell that is frontal in
one step out of twenty happened to catch a passing filament; a cell that
is frontal in fifteen sits on a persistent feature — a shelf-break
front, a tidal mixing front, the edge of a plume — and those are the
ones that aggregate plankton reliably enough for a predator to learn. An
instantaneous distance cannot tell the two apart, and averaging distance
over time does not either, because a cell can be near a different
transient front every step.

## Choosing the threshold

Inherited from
[`distance_to_front()`](https://camilleross.org/derivoce/reference/distance_to_front.md),
and the choice matters more here. With `scope = "record"` one cutoff
applies throughout, so frequency reflects both how often a front is
present and whether this part of the domain is gradient- rich at all.
With `scope = "step"` each step is cut at its own quantile, so a fixed
fraction of cells is frontal in every step and frequency becomes purely
about location. The second is usually what "persistence" is meant to
mean.

## The window

`n = NULL`, the default, uses the whole record and gives a static map:
one value per cell, repeated on every step. That is a description of the
domain rather than a covariate that varies with the observation.

Giving `n` makes it a trailing window in the manner of
[`rolling_covariate()`](https://camilleross.org/derivoce/reference/rolling_covariate.md),
so it varies through the record and can enter a model alongside
conditions at the time. Steps where the cell's gradient is undefined —
the outermost ring, or a missing value — are left out of the denominator
rather than counted as not frontal.

## See also

[`distance_to_front()`](https://camilleross.org/derivoce/reference/distance_to_front.md),
[`rolling_covariate()`](https://camilleross.org/derivoce/reference/rolling_covariate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# A static map of where thermal fronts persist.
env <- front_frequency(env, "SST", scope = "step")

# How frontal this place has been over the preceding year.
env <- front_frequency(env, "SST", scope = "step", n = 12, by = "month")
} # }
```
