# Current speed from velocity components

The magnitude of the horizontal velocity vector, `sqrt(u^2 + v^2)`.
Reproduces the `uv` covariate of Ross et al. (2023).

## Usage

``` r
current_speed(env_dat, u = "UO", v = "VO", name = "speed")
```

## Arguments

- env_dat:

  an `sf` POINT object with one row per location and time step, as
  datamatch's access functions return

- u:

  name of the eastward velocity column

- v:

  name of the northward velocity column

- name:

  name for the new column

## Value

`env_dat` with a current speed column added, in the units of `u`/`v`

## Details

Their `uv_grad` is the spatial derivative of *that speed field*, not of
the velocity components, so it is two steps:

    env <- current_speed(env)                    # uv
    env <- horizontal_gradient(env, "speed")     # uv_grad

The order matters. Differentiating `u` and `v` separately and combining
afterwards gives a different quantity: speed is a non-linear function of
the components, so the gradient of the speed is not the speed of the
gradients.

## References

Ross C, Runge J, Roberts J, Brady D, Tupper B, Record N (2023).
Estimating North Atlantic right whale prey based on Calanus finmarchicus
thresholds. *Marine Ecology Progress Series* **703**, 1-16.
[doi:10.3354/meps14204](https://doi.org/10.3354/meps14204)

## Examples

``` r
if (FALSE) { # \dontrun{
# UO and VO are the defaults, so a dictionary fetch needs no column names
env <- current_speed(env) |> horizontal_gradient("speed")

# If the fetch used Copernicus codes instead
env <- current_speed(env, u = "uo", v = "vo")
} # }
```
