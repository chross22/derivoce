# Describe derived columns as EML metadata

Builds the attribute table that Ecological Metadata Language wants, for
the columns this package added. Derived covariates are the hardest part
of a dataset to document, because their meaning lives in how they were
computed rather than in what was measured: `SST_grad` is degrees per
kilometre by central differences on a lon/lat lattice, and nothing about
the column name or its numbers says so. That knowledge is already in
this package, so it may as well be emitted in the form an archive can
read.

## Usage

``` r
eml_attributes(env_dat, vars = NULL, units = NULL)
```

## Arguments

- env_dat:

  an `sf` POINT object that has been through this package

- vars:

  columns to describe, or `NULL` for every covariate column

- units:

  named character vector giving the EML unit of source columns, for
  example `c(TEMP_INSITU = "celsius")`. Merged over the defaults for the
  variables datamatch serves, which are already known

## Value

a data frame with one row per column, carrying `attributeName`,
`attributeDefinition`, `measurementScale`, `domain`, `unit` and
`numberType`

## Details

The result is a data frame in the shape `EML::set_attributes()`
consumes. derivoce does not depend on the `EML` package and does not
write XML — it hands over the table, and the EML package turns it into a
document:

    attributes <- eml_attributes(env)
    EML::set_attributes(attributes, col_classes = eml_col_classes(env))

## Units it cannot know

A derived unit is usually a function of the source unit — a gradient of
temperature is degrees per kilometre, a gradient of chlorophyll is
mg/m^3 per kilometre — so the source unit has to be known before the
derived one can be.

The variables `datamatch` serves are known already: the Copernicus
physics and biogeochemistry variables, the seafloor terrain from
`attach_bathymetry()`, and the climate indices from
`attach_climate_index()`. A workflow built on those needs no `units`
argument. Pass one for anything else, or to override a default.

Anything still unresolved comes back with `unit = NA` and an
`attributeDefinition` that names the gap, rather than a plausible guess:
an archived dataset with confidently wrong units is worse than one with
an obvious hole in it.

Source columns that this package did not create are described as
`"not derived by derivoce"` for the same reason. They are yours to
document.

## Custom units

EML validates units against a standard dictionary of 195 entries, and
several of the quantities here are not in it — per second, per second
squared, per day, and metres squared per second squared among them.
Those must be declared alongside the attributes or the document will not
validate.
[`eml_custom_units()`](https://camilleross.org/derivoce/reference/eml_custom_units.md)
returns exactly the declarations the attribute table needs, and returns
none if every unit used happened to be standard.

## See also

[`eml_custom_units()`](https://camilleross.org/derivoce/reference/eml_custom_units.md),
[`derived_indices()`](https://camilleross.org/derivoce/reference/derived_indices.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- horizontal_gradient(env, "SST")
env <- eke(env)

attributes <- eml_attributes(env, units = c(SST = "celsius"))
custom <- eml_custom_units(attributes)
} # }
```
