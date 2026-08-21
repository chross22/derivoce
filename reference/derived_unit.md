# A derived unit built from a source unit

Returns `NA` when the source unit is unknown, rather than guessing. A
unit is the one piece of metadata that is worse to state confidently and
wrongly than to leave blank.

## Usage

``` r
derived_unit(source, units, per = NULL)
```

## Arguments

- source:

  the source column name

- units:

  named vector of source units

- per:

  what the source unit is divided by, as an EML-style suffix

## Value

an EML unit name, or `NA`
