# Find eddies as objects, not just as a field

[`flow_deformation()`](https://camilleross.org/derivoce/reference/flow_deformation.md)
returns the Okubo-Weiss parameter as a number per cell. This takes the
next step: it groups the connected cells where rotation beats strain
into individual eddies, and describes each one. A cell then carries not
"how eddy-like is the flow here" but "you are inside an eddy, it turns
this way, and it is this big".

## Usage

``` r
detect_eddies(
  env_dat,
  u = "UO",
  v = "VO",
  threshold = 0.2,
  min_cells = 4L,
  measures = c("in_eddy", "polarity", "radius"),
  per = c("km", "m"),
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

- threshold:

  multiple of the standard deviation of Okubo-Weiss below which a cell
  counts as rotational. Positive; 0.2 follows the literature

- min_cells:

  smallest connected region kept, in cells

- measures:

  any of `"in_eddy"`, `"polarity"` and `"radius"`

- per:

  distance unit for the radius, `"km"` or `"m"`

- suffix:

  appended to each measure to name its column

## Value

`env_dat` with one column per requested measure. `in_eddy` is 1, 0 or
`NA`; `polarity` and `radius` are `NA` outside an eddy

## Details

That distinction matters ecologically, because the two polarities do
opposite things. A cyclonic eddy upwells at its core, lifting nutrients
and often concentrating plankton; an anticyclonic one downwells, and its
core is typically poorer. A covariate that only says "eddy" averages the
two together and can easily find nothing.

## How eddies are identified

The Okubo-Weiss criterion of Isern-Fontanet et al. (2003): a cell
belongs to an eddy where \\W \< -\alpha\sigma_W\\, with \\\sigma_W\\ the
standard deviation of \\W\\ over that time step and \\\alpha\\ the
`threshold` argument. The published value is 0.2, which is the default.
It is a *relative* threshold, recomputed per step, so a quiet month
still yields eddies — the criterion asks which parts of this flow are
most rotational, not whether the flow is energetic in absolute terms.

Connected regions are then labelled, and those smaller than `min_cells`
are discarded as noise. On a 1/12-degree product an eddy of
oceanographic interest spans many cells; a two-cell patch is more likely
to be a numerical artefact of the differencing than a feature.

## What each measure is

- **in_eddy** is 1 inside an identified eddy and 0 outside.

- **polarity** is +1 for cyclonic rotation and -1 for anticyclonic, from
  the sign of the eddy's mean vorticity. In the northern hemisphere
  cyclonic is counter-clockwise and upwelling at the core.

- **radius** is the equivalent radius of the eddy, \\\sqrt{A/\pi}\\ for
  its area \\A\\, in `per` units. Real eddies are not discs, so this is
  a size, not a shape. See below for what it does and does not measure.

Cells outside any eddy get 0 for `in_eddy` and `NA` for the others,
because they have no eddy to describe. Cells where there was nothing to
judge get `NA` for all three: the outermost ring, where the central
difference has no neighbours, and anywhere the velocity itself is
missing, which over a masked product is land. Those are not the same as
being outside an eddy, and `in_eddy` does not report them as 0 – a
fabricated zero over land is indistinguishable from a measured one, and
would be used as though it were.

## Radius measures the patch, not the eddy

The threshold is a multiple of \\\sigma_W\\ for the whole time step, so
it is one absolute level applied to every eddy in the field at once. A
strong eddy has \\\|W\|\\ far below that level over nearly its whole
core and loses almost nothing to it; a weak one only dips below it near
its centre, and the patch that survives is a fraction of the feature.
Two eddies of the same physical size can therefore report very different
radii if one is spinning faster than the other.

On a field of Gaussian vortices with cores of 25 to 60 km, raising
`threshold` from 0.05 to 1.5 moves the reported radius of the strongest
eddy by 3 percent and that of the weakest by more than half. `radius` is
comparable between eddies of similar strength, and between time steps
only where \\\sigma_W\\ is stable. Where the question is size as such,
treat it as an index rather than a measurement, and hold `threshold`
fixed across everything being compared.

## What this is not

Detection per time step, with no identity between steps. Nothing here
tracks an eddy through its life, so there is no age, no lifespan and no
propagation speed, and the same physical eddy carries unrelated labels
in consecutive steps. `eddy_id` is therefore deliberately not returned:
it would invite exactly that mistake.

The Okubo-Weiss criterion is also known to be permissive in strongly
strained flow and sensitive to the smoothness of the velocity field. It
finds rotation-dominated regions, which is a good working definition of
an eddy and not the same thing as a closed-contour eddy from altimetry.

## References

Isern-Fontanet J, Garcia-Ladona E, Font J (2003). Identification of
marine eddies from altimetric maps. *Journal of Atmospheric and Oceanic
Technology* **20**(5), 772-778.
[doi:10.1175/1520-0426(2003)20\<772:IOMEFA\>2.0.CO;2](https://doi.org/10.1175/1520-0426%282003%2920%3C772%3AIOMEFA%3E2.0.CO%3B2)

## See also

[`flow_deformation()`](https://camilleross.org/derivoce/reference/flow_deformation.md),
[`distance_to_eddy()`](https://camilleross.org/derivoce/reference/distance_to_eddy.md),
[`eke()`](https://camilleross.org/derivoce/reference/eke.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- detect_eddies(env)

# Cyclonic cores upwell, so they are the ones worth separating out.
upwelling <- env$in_eddy == 1 & env$polarity > 0
} # }
```
