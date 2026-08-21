# Regional indices this package can derive

The named region-scale indices, what each measures, what it needs as
input, and where the concept comes from. All are derived from gridded
fields rather than downloaded, which distinguishes them from the climate
indices
[`datamatch::attach_climate_index()`](https://camilleross.org/datamatch/reference/attach_climate_index.html)
serves.

## Usage

``` r
derived_indices(markdown = FALSE)

# S3 method for class 'derivoce_index_table'
print(x, ...)
```

## Arguments

- markdown:

  return a markdown table instead of printing a formatted one, for
  pasting into documentation

- x:

  a `derivoce_index_table`

- ...:

  ignored

## Value

a data frame of class `derivoce_index_table`, or a character string of
markdown when `markdown = TRUE`

## Why three ways of measuring the same thing

Scotian Shelf inflow to the Gulf of Maine appears here three times,
because the literature measures it three ways and they answer different
questions.

A **transport** across a section measures the crossing itself, and is
the only one that gives a direction and a flux. A **water-mass
fraction** measures how much of the water present came from there, which
is what matters for nutrients and is available even where velocities are
not. A **box anomaly** measures the most visible consequence,
freshening, and is the most robust to poor velocity data but the least
specific about cause.

They disagree in informative ways. A strong inflow with a normal
salinity anomaly means the arriving water was not unusually fresh, which
is a real finding about the upstream shelf rather than a contradiction.

## See also

[`scotian_shelf_inflow()`](https://camilleross.org/derivoce/reference/scotian_shelf_inflow.md),
[`northeast_channel_inflow()`](https://camilleross.org/derivoce/reference/northeast_channel_inflow.md),
[`water_mass_fraction()`](https://camilleross.org/derivoce/reference/water_mass_fraction.md),
[`eastern_gom_salinity()`](https://camilleross.org/derivoce/reference/eastern_gom_salinity.md)

## Examples

``` r
derived_indices()
#> Regional indices derivoce can compute
#> ----------------------------------------------------------------------
#>  name                     method               needs                       
#>  scotian_shelf_inflow     transport            UO, VO                      
#>  northeast_channel_inflow transport            UO, VO                      
#>  water_mass_fraction      T-S endmember mixing SST, SSS, and two endmembers
#>  eastern_gom_salinity     box anomaly          SSS                         
#>  units           
#>  m^2/s           
#>  m^2/s           
#>  fraction, 0 to 1
#>  PSU             
#> 
#> scotian_shelf_inflow
#>   measures: Scotian Shelf Water crossing into the Gulf of Maine at Cape Sable
#>   sign:     positive into the Gulf of Maine
#>   source:   Feng et al. 2016; Wang et al. 2022
#> 
#> northeast_channel_inflow
#>   measures: slope water entering the Gulf of Maine through the Northeast Channel
#>   sign:     positive into the Gulf of Maine
#>   source:   Ramp et al. 1985; Du et al. 2022; Silver et al. 2023
#> 
#> water_mass_fraction
#>   measures: what fraction of each cell is a named water mass
#>   sign:     1 is entirely the first endmember
#>   source:   Townsend et al. 2015
#> 
#> eastern_gom_salinity
#>   measures: freshening in the eastern Gulf of Maine, where inflow first appears
#>   sign:     negative is fresher, indicating stronger inflow
#>   source:   Grodsky et al. 2025
#> 
#> Full descriptions: as.data.frame(derived_indices())$description
#> Markdown table:    cat(derived_indices(markdown = TRUE))

cat(derived_indices(markdown = TRUE))
#> | Index | Measures | Method | Needs | Units | Source |
#> | --- | --- | --- | --- | --- | --- |
#> | `scotian_shelf_inflow()` | Scotian Shelf Water crossing into the Gulf of Maine at Cape Sable | transport | UO, VO | m^2/s | Feng et al. 2016; Wang et al. 2022 |
#> | `northeast_channel_inflow()` | slope water entering the Gulf of Maine through the Northeast Channel | transport | UO, VO | m^2/s | Ramp et al. 1985; Du et al. 2022; Silver et al. 2023 |
#> | `water_mass_fraction()` | what fraction of each cell is a named water mass | T-S endmember mixing | SST, SSS, and two endmembers | fraction, 0 to 1 | Townsend et al. 2015 |
#> | `eastern_gom_salinity()` | freshening in the eastern Gulf of Maine, where inflow first appears | box anomaly | SSS | PSU | Grodsky et al. 2025 |
```
