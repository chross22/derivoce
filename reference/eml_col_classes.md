# Column classes for EML, in the order the attributes are given

`EML::set_attributes()` wants a `col_classes` vector alongside the
table, one entry per attribute, drawn from "numeric", "character",
"factor" and "Date".

## Usage

``` r
eml_col_classes(env_dat, vars = NULL)
```

## Arguments

- env_dat:

  an `sf` POINT object

- vars:

  columns, or `NULL` for every covariate column

## Value

a character vector

## See also

[`eml_attributes()`](https://camilleross.org/derivoce/reference/eml_attributes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
EML::set_attributes(eml_attributes(env), col_classes = eml_col_classes(env))
} # }
```
