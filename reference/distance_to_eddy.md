# Distance to the nearest eddy

The sibling of
[`distance_to_front()`](https://camilleross.org/derivoce/reference/distance_to_front.md)
and
[`distance_to_isobath()`](https://camilleross.org/derivoce/reference/distance_to_isobath.md).
Being outside an eddy but close to one is a different place from being
far from any, and a binary inside/outside flag cannot express it: a
station 5 km from a rotating core and one 300 km away both score zero.

## Usage

``` r
distance_to_eddy(
  env_dat,
  u = "UO",
  v = "VO",
  threshold = 0.2,
  min_cells = 4L,
  polarity = c("any", "cyclonic", "anticyclonic"),
  per = c("km", "m"),
  name = NULL
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

- polarity:

  `"any"`, `"cyclonic"`, or `"anticyclonic"`

- per:

  distance unit for the radius, `"km"` or `"m"`

- name:

  name for the new column

## Value

`env_dat` with a distance column, in `per` units. A cell inside a
qualifying eddy has a distance of 0. A step containing no qualifying
eddy returns `NA` rather than a distance to nothing

## Details

`polarity` narrows it to one kind. Because cyclonic and anticyclonic
eddies do opposite things to the water column, distance to the nearest
*cyclonic* eddy is often the covariate that carries a signal where
distance to the nearest eddy of any sort does not.

## See also

[`detect_eddies()`](https://camilleross.org/derivoce/reference/detect_eddies.md),
[`distance_to_front()`](https://camilleross.org/derivoce/reference/distance_to_front.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- distance_to_eddy(env, polarity = "cyclonic")
} # }
```
