# What each suffixed derived column means

Matched against the end of a column name, longest patterns first so that
`_front_dist` is not mistaken for a rolling statistic. Each entry
supplies a function giving the unit, because most derived units are
built from the source column's own unit.

## Usage

``` r
eml_suffix_registry()
```

## Value

a list of pattern entries
