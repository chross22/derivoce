# Every column a deposited table should describe

Wider than
[`covariate_columns()`](https://camilleross.org/derivoce/reference/covariate_columns.md),
and deliberately so. `<var>_depth` must not be derived from or resampled
as though it were a measurement, but an archive still has to say what it
is -
[`datamatch::write_eml()`](https://camilleross.org/datamatch/reference/write_eml.html)
documents it, with metres as its unit, rather than leaving a column in
the deposited table with nothing describing it. Only datamatch's
internal per-row tag is left out, since it does not survive
`matchData()` to reach a deposit at all.

## Usage

``` r
documented_columns(env_dat)
```

## Arguments

- env_dat:

  an `sf` POINT object

## Value

character vector of column names
