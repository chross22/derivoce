# Distance to the nearest front

Fronts are where a covariate changes sharply over a short distance, and
they concentrate plankton. But the local gradient alone is a blunt
predictor of that: a station sitting in the middle of a smooth patch has
a gradient of zero whether the nearest front is 2 km away or 200 km
away, and those are very different places to be. Distance to the nearest
front separates them.

## Usage

``` r
distance_to_front(
  env_dat,
  var,
  threshold = NULL,
  quantile = 0.9,
  scope = c("record", "step"),
  per = c("km", "m"),
  name = NULL
)
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- var:

  covariate whose gradient defines the fronts

- threshold:

  absolute gradient value above which a cell is frontal; takes
  precedence over `quantile` when given

- quantile:

  fraction of the gradient distribution treated as frontal

- scope:

  `"record"` or `"step"`; see above

- per:

  distance unit for the result: `"km"` (default) or `"m"`

- name:

  name for the new column

## Value

`env_dat` with a distance-to-front column. Cells on a front are `0`. A
time step containing no front is `NA` throughout, since "distance to
nothing" has no value.

## Details

Fronts are identified by thresholding the covariate's horizontal
gradient, then the distance from every cell to the nearest front cell is
computed with a distance transform.

## Choosing the threshold

`quantile` (the default route) calls a cell frontal if its gradient is
in the top fraction of the distribution — `0.9` means the strongest 10%
of gradients.

`scope` decides what that fraction is taken over, and the two answers
mean different things:

- `"record"` (the default) takes one threshold across the whole time
  series, so "front" means the same physical sharpness in every month.
  Winter months with weak stratification may then contain no fronts at
  all — which is a real result, not a failure.

- `"step"` takes a separate threshold per time step, so every month has
  fronts by construction. Use this when the question is where the
  sharpest features are *this month*, and not whether this month has
  sharp features.

Pass `threshold` instead to set an absolute gradient value, which is the
right choice when a physically meaningful cutoff is known.

## References

Identifying fronts by thresholding a gradient field follows the standard
approach, of which Belkin and O'Reilly's is the reference implementation
for satellite SST and chlorophyll. Theirs adds a contextual median
filter that preserves front shape before the gradient is taken, which
matters on noisy satellite retrievals and less on a model field that is
already smooth. What is done here is the plain gradient threshold, with
`scope` deciding whether the threshold is one value for the record or
one per time step.

Belkin IM, O'Reilly JE (2009). An algorithm for oceanic front detection
in chlorophyll and SST satellite imagery. *Journal of Marine Systems*
**78**(3), 319-326.
[doi:10.1016/j.jmarsys.2008.11.018](https://doi.org/10.1016/j.jmarsys.2008.11.018)

## Examples

``` r
if (FALSE) { # \dontrun{
# Distance to the sharpest 10% of thermal gradients
env <- distance_to_front(env, "SST")
# A physically defined front: 0.05 degrees C per km
env <- distance_to_front(env, "SST", threshold = 0.05)
} # }
```
