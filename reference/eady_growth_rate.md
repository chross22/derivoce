# Eady growth rate of baroclinic instability

How fast baroclinic instability grows: the rate at which the available
potential energy stored in tilted density surfaces is converted into
eddies. High where a sheared, weakly stratified flow can overturn, which
on this shelf means the shelf-break front and the edges of warm-core
rings — the places eddies are actually generated.

## Usage

``` r
eady_growth_rate(
  env_dat,
  shallow = c("UO", "VO"),
  deep,
  depths,
  stratification = "N2",
  per = c("day", "second"),
  name = NULL
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- shallow:

  length-2 character vector naming the eastward and northward velocity
  columns at the shallower level, in m/s

- deep:

  the same for the deeper level

- depths:

  the two depths in metres, shallower first

- stratification:

  name of an \\N^2\\ column, in s^-2, as produced by
  [`buoyancy_frequency()`](https://camilleross.org/derivoce/reference/buoyancy_frequency.md)

- per:

  `"day"` (the default) or `"second"`

- name:

  name for the new column

## Value

`env_dat` with a growth rate column, in day^-1 or s^-1

## Details

Computed in the maximum-growth-rate form of Lindzen and Farrell (1980),
\\\sigma = 0.31\\\|f\|\\\|\partial U/\partial z\| / N\\, for the
instability described by Eady (1949).

**Eady is a person, not a spelling of "eddy".** Eric Eady set out the
model in 1949. The collision is unlucky, because the rate named after
him is precisely a predictor of where eddies form, so the two words look
interchangeable beside each other and are not.

## What it needs

Velocities at two depths and a stratification spanning the same layer.
From Copernicus that means two calls at different `depth` ranges joined
as columns; from HYCOM it means `UO`/`VO` beside
`UO_BOTTOM`/`VO_BOTTOM`, which arrive together. See
[`buoyancy_frequency()`](https://camilleross.org/derivoce/reference/buoyancy_frequency.md),
which produces the stratification and explains how to choose the levels.

## What it is and is not

An index of **where and when** instability is favoured, not a prediction
of growth. The formula assumes quasi-geostrophic scaling and uniform
stratification through the layer, neither of which holds exactly on a
shelf, and the 0.31 is the maximum over wavenumber rather than the rate
of any particular disturbance. Read it as a comparative field.

It is undefined where the column is not stably stratified, since \\N\\
is then not a frequency. Those cells come back `NA` rather than as a
very large rate, which is what dividing by a vanishing \\N\\ would
otherwise produce and which would look like intense instability exactly
where the assumption has failed.

## References

Eady ET (1949). Long waves and cyclone waves. *Tellus* **1**(3), 33-52.
[doi:10.3402/tellusa.v1i3.8507](https://doi.org/10.3402/tellusa.v1i3.8507)

Lindzen RS, Farrell B (1980). A simple approximate result for the
maximum growth rate of baroclinic instabilities. *Journal of the
Atmospheric Sciences* **37**(7), 1648-1654.
[doi:10.1175/1520-0469(1980)037\<1648:ASARFT\>2.0.CO;2](https://doi.org/10.1175/1520-0469%281980%29037%3C1648%3AASARFT%3E2.0.CO%3B2)

## See also

[`buoyancy_frequency()`](https://camilleross.org/derivoce/reference/buoyancy_frequency.md),
[`detect_eddies()`](https://camilleross.org/derivoce/reference/detect_eddies.md),
[`flow_deformation()`](https://camilleross.org/derivoce/reference/flow_deformation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- buoyancy_frequency(env, "rho_surface", "rho_deep",
                          depths = c(0.494, 92.326))
env <- eady_growth_rate(env,
                        shallow = c("UO", "VO"),
                        deep = c("UO_deep", "VO_deep"),
                        depths = c(0.494, 92.326))
} # }
```
