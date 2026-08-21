# Is a column constant within every group?

The shared test behind both degeneracy warnings: group by location and a
constant column is static through time, group by time step and it is
uniform across space.

## Usage

``` r
constant_within(values, group)
```

## Arguments

- values:

  a covariate column

- group:

  grouping vector of the same length

## Value

`TRUE` if every group holds a single distinct non-`NA` value

## Details

Implemented by counting distinct (group, value) pairs rather than by
splitting and looping, so it stays O(n) on the full point table.
Comparison is exact: two values differing in the last bit count as
different, so this only fires on columns that really are constant, never
on a nearly-flat field.
