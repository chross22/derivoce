# Velocity gradient diagnostics: vorticity, divergence, strain, Okubo-Weiss

Derives the local structure of a flow from the four horizontal
derivatives of its velocity components. Together these say whether a
parcel of water is spinning, converging, being pulled apart, or some
combination.

## Usage

``` r
flow_deformation(
  env_dat,
  u = "UO",
  v = "VO",
  measures = c("vorticity", "divergence", "strain_rate", "okubo_weiss"),
  suffix = ""
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- u:

  name of the eastward velocity column, in m/s

- v:

  name of the northward velocity column, in m/s

- measures:

  which diagnostics to add. Any of `"vorticity"`, `"divergence"`,
  `"normal_strain"`, `"shear_strain"`, `"strain_rate"`, `"okubo_weiss"`
  and `"rossby"`

- suffix:

  appended to each measure to name its column

## Value

`env_dat` with one column per requested measure. Vorticity, divergence
and the strains are in s^-1, Okubo-Weiss in s^-2, and the Rossby number
is dimensionless

## Details

All are instantaneous and local: they describe the flow in one cell at
one time step, and need no history. That places them between
[`eke()`](https://camilleross.org/derivoce/reference/eke.md), which
needs a series to form an anomaly, and
[`ftle()`](https://camilleross.org/derivoce/reference/ftle.md) and
[`fsle()`](https://camilleross.org/derivoce/reference/fsle.md), which
advect particles through many steps. Where those two answer "how
variable is this place" and "where do trajectories separate", these
answer "what is the flow doing here, right now".

## The measures

Writing the derivatives of eastward `u` and northward `v` with respect
to eastward `x` and northward `y`:

- **vorticity**, \\\zeta = \partial v/\partial x - \partial u/\partial
  y\\, is rotation. Positive is counter-clockwise (cyclonic in the
  northern hemisphere), which in an eddy means upwelling at its core.

- **divergence**, \\\partial u/\partial x + \partial v/\partial y\\, is
  spreading. Negative is convergence, where surface water piles up and
  sinks, and where buoyant particles and floating material accumulate.

- **normal_strain**, \\\partial u/\partial x - \partial v/\partial y\\,
  and **shear_strain**, \\\partial v/\partial x + \partial u/\partial
  y\\, stretch a parcel along the axes and diagonally.

- **strain_rate** is their magnitude, the total rate of deformation.

- **okubo_weiss**, \\W = S_n^2 + S_s^2 - \zeta^2\\, compares the two.
  Negative means rotation wins, which is the interior of a coherent
  eddy. Positive means strain wins, which is the filaments between
  eddies where water is drawn out into long thin structures.

- **rossby** is \\\zeta/f\\, vorticity as a fraction of planetary
  rotation. It is the dimensionless version, so it can be compared
  between latitudes, and values approaching 1 mean the flow is fast
  enough for the usual geostrophic balance to be breaking down.

## What it cannot tell you

These are single-time-step quantities, so a large value means the flow
is deforming now, not that it has been or will be. A coherent eddy that
persists for weeks and a momentary filament can carry the same
Okubo-Weiss value. Persistence is what
[`ftle()`](https://camilleross.org/derivoce/reference/ftle.md) and
[`fsle()`](https://camilleross.org/derivoce/reference/fsle.md) measure,
at much greater cost.

The derivatives are central differences, so like every gradient in this
package the outermost ring of cells is `NA`, and the values are only as
smooth as the grid. On a coarse product a real eddy smaller than a few
cells will not appear at all.

## References

Okubo A (1970). Horizontal dispersion of floatable particles in the
vicinity of velocity singularities such as convergences. *Deep-Sea
Research and Oceanographic Abstracts* **17**(3), 445-454.
[doi:10.1016/0011-7471(70)90059-8](https://doi.org/10.1016/0011-7471%2870%2990059-8)

Weiss J (1991). The dynamics of enstrophy transfer in two-dimensional
hydrodynamics. *Physica D: Nonlinear Phenomena* **48**(2-3), 273-294.
[doi:10.1016/0167-2789(91)90088-Q](https://doi.org/10.1016/0167-2789%2891%2990088-Q)

Isern-Fontanet J, Garcia-Ladona E, Font J (2003). Identification of
marine eddies from altimetric maps. *Journal of Atmospheric and Oceanic
Technology* **20**(5), 772-778.
[doi:10.1175/1520-0426(2003)20\<772:IOMEFA\>2.0.CO;2](https://doi.org/10.1175/1520-0426%282003%2920%3C772%3AIOMEFA%3E2.0.CO%3B2)

## See also

[`eke()`](https://camilleross.org/derivoce/reference/eke.md),
[`ftle()`](https://camilleross.org/derivoce/reference/ftle.md),
[`fsle()`](https://camilleross.org/derivoce/reference/fsle.md),
[`horizontal_gradient()`](https://camilleross.org/derivoce/reference/horizontal_gradient.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- flow_deformation(env)

# Eddy interiors are where rotation beats strain.
eddy_core <- env$okubo_weiss < 0
} # }
```
