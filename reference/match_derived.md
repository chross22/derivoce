# Match a column name against what this package produces

Derived columns are named by convention rather than registered, so the
name is what identifies them. Fixed names are matched outright and
suffixes are matched against the stem, which also yields the source
column whose unit the derived unit is built from.

## Usage

``` r
match_derived(name, units)
```

## Arguments

- name:

  the column name

- units:

  named vector of source units

## Value

a list with `definition`, `unit` and `number`
