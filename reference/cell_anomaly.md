# Anomaly of a covariate against each cell's own history

Subtracts, from every value, the mean of that same grid cell. The result
is a **map per time step** saying where conditions departed from what
that place is usually like, rather than one number per step.

## Usage

``` r
cell_anomaly(
  env_dat,
  vars = NULL,
  reference = c("climatology", "record"),
  standardize = FALSE,
  detrend = FALSE,
  suffix = NULL
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- vars:

  covariate columns, or `NULL` for all numeric ones

- reference:

  `"climatology"` (the default) removes a separate mean per cell per
  calendar month, leaving departures from the usual conditions for the
  time of year. `"record"` removes one mean per cell over the whole
  series, which leaves the seasonal cycle in

- standardize:

  divide by the standard deviation of the same group, giving a z-score

- detrend:

  also remove a fitted linear trend from each cell, so what is left is
  the departure after both the seasonal cycle and the long-term change
  are gone

- suffix:

  appended to each covariate name to make the new column. Defaults to
  `"_anom"`, or `"_z"` when standardising

## Value

`env_dat` with one anomaly column per covariate, in the units of the
covariate, or dimensionless when standardised

## Details

This is the cell-wise counterpart of
[`box_anomaly()`](https://camilleross.org/derivoce/reference/box_anomaly.md).
Where that averages a region into an index, this leaves the spatial
pattern intact and removes the spatial pattern of the *mean* instead. It
is often what makes a covariate usable across a domain with a strong
background gradient: 8 degrees is cold for the southern Gulf of Maine
and warm for the Scotian Shelf, and a model given raw temperature has to
learn that geography before it can use the departure. An anomaly hands
it over directly.

## Standardising

`standardize = TRUE` divides by the cell's standard deviation as well,
giving a z-score. That makes departures comparable between places with
different variability — a 1 degree anomaly is unremarkable on the shelf
and extreme in the deep basin — which is what you want when a single
coefficient has to apply across the whole domain.

It also throws away the magnitude. If the question is how much warmer,
not how unusual, leave it off.

## Detrending

`detrend = TRUE` removes a fitted linear trend as well. This matters in
a warming shelf sea: an anomaly that still contains the trend largely
encodes *which year it is*, and a model given it will fit the trend and
appear to have learned something about temperature. What is left after
detrending says whether conditions were unusual **for their year**,
which is usually the ecological question.

The trend and the seasonal cycle are estimated in one fit rather than
sequentially, for the reasons set out in
[`decompose_covariate()`](https://camilleross.org/derivoce/reference/decompose_covariate.md)
— removed one after the other, each absorbs part of the other. With
`reference = "climatology"` both are taken out and the result is that
function's residual; with `reference = "record"` only the mean and the
trend go, and the seasonal cycle stays in.

Detrending is not free of assumptions. A linear trend fitted to a short
record can absorb genuine low-frequency variability — a decade of a
multidecadal oscillation looks like a trend — so the residual is a
departure from a fitted line, not from a known baseline.

## Why this needs several years

With `reference = "climatology"` the mean is taken over the values
sharing a calendar month, so a cell's January is compared with its other
Januaries. Given a single year of data there is exactly one January per
cell, the mean is that value, and **every anomaly is identically zero**.
That is a silent failure — a column of zeroes is a perfectly
plausible-looking covariate — so it warns.

The same applies to `standardize`, which needs at least two values per
group to have a standard deviation at all.

## See also

[`box_anomaly()`](https://camilleross.org/derivoce/reference/box_anomaly.md),
[`marine_heatwave()`](https://camilleross.org/derivoce/reference/marine_heatwave.md),
[`lag_covariate()`](https://camilleross.org/derivoce/reference/lag_covariate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Departures from the usual conditions for the month, cell by cell.
env <- cell_anomaly(env, "SST")

# Comparable across a domain with very different variability.
env <- cell_anomaly(env, "SST", standardize = TRUE)
} # }
```
