# Declarations for the units EML does not know

EML validates against a fixed dictionary, and a unit outside it must be
declared as a `customUnit` or the document fails validation. Several
quantities in this package are outside it, so this returns the
declarations for whichever ones an attribute table actually used.

## Usage

``` r
eml_custom_units(attributes)
```

## Arguments

- attributes:

  a table from
  [`eml_attributes()`](https://camilleross.org/derivoce/reference/eml_attributes.md),
  or a character vector of unit names

## Value

a data frame of custom unit declarations, with no rows if every unit
used is standard

## Details

Each row carries what EML requires of a custom unit: an `id` matching
the `unit` in the attribute table, the `unitType` it belongs to, the SI
unit it derives from, the multiplier to reach SI, and a description.

## See also

[`eml_attributes()`](https://camilleross.org/derivoce/reference/eml_attributes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
attributes <- eml_attributes(env)
EML::set_unitList(eml_custom_units(attributes))
} # }
```
