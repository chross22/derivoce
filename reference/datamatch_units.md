# EML units for the variables datamatch serves

A translation of datamatch's own unit strings into EML unit ids, so that
[`eml_attributes()`](https://camilleross.org/derivoce/reference/eml_attributes.md)
can resolve the source columns of a normal workflow without being told
what each one holds. Derived units are built on top of these, so a
temperature gradient becomes degrees per kilometre without the caller
naming either part.

## Usage

``` r
datamatch_units()
```

## Value

a named character vector, variable name to EML unit id

## Details

Hardcoded rather than read from datamatch at runtime, because datamatch
is deliberately not a dependency of this package. The cost is that it
can fall behind, so a test checks it against the live catalogue whenever
datamatch is installed, and reports any variable this table has not been
told about.
