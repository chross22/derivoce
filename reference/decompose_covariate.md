# Split a covariate into trend, seasonal cycle, and what is left

Separates each cell's series into the parts it is made of: a long-term
trend, a repeating seasonal cycle, and the residual departure from both.
The pieces are additive, so `value = mean + trend + seasonal + residual`
to numerical precision, and any of them can be used on its own.

## Usage

``` r
decompose_covariate(
  env_dat,
  vars = NULL,
  degree = 1,
  components = c("trend", "seasonal", "residual", "slope"),
  suffix = "_"
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- vars:

  covariate columns, or `NULL` for all numeric ones

- degree:

  polynomial degree of the trend. 1 is a straight line

- components:

  which of `"trend"`, `"seasonal"`, `"residual"` and `"slope"` to add

- suffix:

  appended to each covariate name, before the component name

## Value

`env_dat` with one column per covariate per component

## Details

The separation is worth making because the parts answer different
questions and are easy to confuse. In a warming shelf sea, a raw
temperature anomaly that still contains the trend largely encodes *which
year it is*: a model given it will fit the trend and appear to have
learned something about temperature. The residual is the part that says
whether this month was warm *for its year and season*, which is usually
the ecological question.

## What each part is

- **trend** is the fitted long-term component, centred on zero, in the
  units of the covariate. `degree = 1` is a straight line; raise it
  where the trend curves, and be aware that a high degree will start
  absorbing variability that is not a trend.

- **seasonal** is the average departure for each calendar month, taken
  after the trend is removed and centred on zero. It repeats every year
  by construction, so it cannot represent a seasonal cycle that is
  itself changing.

- **residual** is everything else, which is what
  [`cell_anomaly()`](https://camilleross.org/derivoce/reference/cell_anomaly.md)
  with `detrend = TRUE` returns.

- **slope** is the rate of change per year, one number per cell,
  repeated on its rows. Only defined for `degree = 1`, where the trend
  is a line and a single rate describes it.

## What it needs, and what it does when it does not have it

Fitting a trend of `degree` needs at least `degree + 1` time steps for
that cell, and a monthly climatology needs several years before a
month's mean is anything but that month's one value. A cell with too
little history gets `NA` for the parts that cannot be estimated rather
than a fit through two points, and the function warns.

The trend and the seasonal cycle are estimated **together**, in one fit,
rather than one being removed before the other is measured. Doing it in
sequence lets whichever goes first absorb part of the other, in both
directions.

A trend removed first picks up part of the seasonal cycle, because the
trend is fitted against elapsed days and calendar months are of unequal
length: a twelve-month cycle sampled on the first of each month is not
orthogonal to a straight line in days, and a series containing nothing
but a seasonal cycle comes out with a spurious trend worth a few percent
of its amplitude. A seasonal cycle removed first picks up part of the
trend whenever the record does not contain whole years — with a series
starting in July and warming throughout, the later months carry more
warm years than the earlier ones, and the "cycle" acquires a step that
is really the trend.

Fitting both at once removes each conditional on the other, so neither
failure occurs. The cost is that the seasonal term needs enough
observations to afford a parameter per calendar month; where it cannot,
it is reported as zero rather than estimated badly, and the function
warns.

## See also

[`cell_anomaly()`](https://camilleross.org/derivoce/reference/cell_anomaly.md),
[`box_anomaly()`](https://camilleross.org/derivoce/reference/box_anomaly.md),
[`index_series()`](https://camilleross.org/derivoce/reference/index_series.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- decompose_covariate(env, "SST")
# adds SST_trend, SST_seasonal, SST_residual

# How fast is each cell warming, in degrees per year?
env <- decompose_covariate(env, "SST", components = "slope")
} # }
```
