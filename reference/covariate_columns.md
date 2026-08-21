# Covariate column names in an environmental data object

Deliberately internal, and a deliberate duplicate of
[`datamatch::covariate_columns()`](https://camilleross.org/datamatch/reference/covariate_columns.html),
which is identical. Exporting it would mask datamatch's whenever both
packages are attached, for no gain: callers who want it have it from
datamatch already. It is not imported from there because datamatch is a
suggested rather than a hard dependency — everything in this package
operates on the *shape* datamatch's access functions return, not on
datamatch itself, so an object of that shape from any source works. That
shape is common to all of them – `accessCopernicus()`, `accessHYCOM()`,
`accessCCMP()`, `accessFVCOM()` and `accessERDDAP()` – and to anything
else that produces one row per location and time step.

## Usage

``` r
covariate_columns(env_dat)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

## Value

character vector of covariate column names

## Details

`<var>_source` is included and `<var>_depth` is not, matching
[`datamatch::covariate_columns()`](https://camilleross.org/datamatch/reference/covariate_columns.html).
A source tag travels with the variable it describes and survives
resampling as a categorical; a depth is the model level a derived bottom
value was taken from, and the mean of two of those is not the depth any
value came from. `.datamatch_source` is datamatch's internal per-row tag
and is never data.
