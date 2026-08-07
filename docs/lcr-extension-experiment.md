# Can the Labrador Current retroflection index be extended from Copernicus?

**Short answer: no, on two separate attempts.** Monthly fields fail for a
physical reason. Daily fields with OceanParcels, which is what the first attempt
concluded was needed, clear that obstacle and still do not reproduce the
published series. This records both, so nobody repeats them.

Read the monthly attempt for why daily fields are necessary, and
[the daily attempt](#the-daily-attempt-oceanparcels) for why they are not
sufficient. The pipeline built for the second attempt is in
[`../lcr-pipeline/`](../lcr-pipeline/) and is kept, working and tested, in case
the missing piece ever becomes available.

## Why it was worth trying

`datamatch::attach_climate_index()` serves `LCR`, the retroflection index of
Jutras et al. (2023). It covers **1993–2014 only** and cannot be attached to
recent observations, which is a real limitation for anything modelling the last
decade.

The index is derived rather than observed, so in principle it can be recomputed.
Jutras et al. seeded virtual particles weekly across a line on the Labrador Shelf
at (53°N, 56.7°W)–(54.3°N, 52.0°W), tracked each for three years through
GLORYS12V1 with [OceanParcels](https://oceanparcels.org/), and took the
difference between the counts reaching hydrographic sections on the Labrador and
Scotian shelves. High retroflection means fewer particles continue southwest.

derivoce already advects particles: `ftle()` and `fsle()` use fourth-order
Runge-Kutta through a time-varying field, with the same m/s to degrees-per-day
conversion. Reusing that machinery to recompute the index looked cheap.

## What was done

Monthly GLORYS surface velocities, 1993–2000, over 40–58°N and 70–45°W, which
spans the release line and the Scotian Shelf. Particles seeded along the Jutras
line, advected for three years on a two-day step, counting those that reached a
Labrador Shelf box and a Scotian Shelf box.

The script is [`lcr-extension-experiment.R`](lcr-extension-experiment.R).

## What happened

Particles reached the Labrador Shelf and **never reached the Scotian Shelf**. Not
some months, none at all. They also survived a median of 112 days of the 1095
asked for, leaving the domain to the east.

Tracing where they went explains it: they travelled south but **never west of
56°W**, and the Scotian Shelf begins around 59°W.

The eight-year mean surface flow says why. There is no coherent southwestward
shelf jet in it:

| Location | u (m/s) | v (m/s) | Direction |
|---|---|---|---|
| Release line (54°W, 53.6°N) | −0.015 | −0.077 | southwest |
| Labrador shelf (55°W, 52°N) | **+0.050** | −0.286 | south**east** |
| Newfoundland (52°W, 49°N) | +0.070 | −0.029 | southeast |
| Grand Banks (50°W, 46°N) | +0.008 | −0.003 | near stagnant |
| Scotian Shelf (62°W, 44°N) | +0.033 | −0.019 | east |

The Labrador Current is a narrow shelf-edge jet, of order 100 km wide. A monthly
mean smears it, and at the surface it competes with wind-driven Ekman drift
pointing offshore. What survives averaging is a southeastward drift into the open
Atlantic rather than a current following the shelf.

The Grand Banks row is the fatal one. That is where retroflection is decided, and
in the monthly mean the flow there is essentially zero. There is nothing to steer
a particle either way, so the quantity the index measures does not exist in this
field.

## Depth does not rescue it

The surface layer is the most Ekman-contaminated, so the same comparison was run
at 40 m:

| Location | surface u | 40 m u |
|---|---|---|
| Labrador shelf | +0.050 | +0.018 |
| Grand Banks | +0.008 | −0.012 |

Better, and still not a shelf-following jet. The eastward component on the
Labrador shelf shrinks but does not reverse.

## What would be needed

**Daily fields**, as Jutras et al. used. The eddies and the jet survive daily
averaging and do not survive monthly averaging.

That is a much larger undertaking than a bigger download. For 1993–2014 over this
domain, daily GLORYS is roughly 4 GB. Worse is the computation: weekly releases
over 22 years is about 1,100 releases, each tracking of order a thousand
particles for 1,095 daily steps, with four velocity evaluations per step. That is
on the order of a billion interpolations. In R, with `terra::extract` per step,
the pilot here ran about 11 ms per step for 33 particles, which extrapolates to
roughly a day of compute, and the daily fields cannot all be held in memory.

OceanParcels exists for exactly this, which is presumably why Jutras et al. used
it.

## Conclusion of the monthly attempt

Recomputing `LCR` from monthly Copernicus fields does not work, and the failure
is physical rather than a matter of tuning. The index depends on a narrow current
and a bifurcation that monthly averaging removes.

Anyone wanting to extend the series past 2014 should plan on daily fields and a
purpose-built Lagrangian framework, not on this package's advection.

# The daily attempt: OceanParcels

That recommendation was then followed. Daily GLORYS, OceanParcels, weekly
releases along the Jutras line, three years of tracking, arrivals counted in
Labrador and Scotian boxes. The pipeline is
[`../lcr-pipeline/`](../lcr-pipeline/).

It works. It does not reproduce the index.

## Depth decides everything, and there is a narrow window

The monthly attempt already suspected Ekman drift. With daily fields the effect
is sharp, and a second effect appears that is easy to miss: **GLORYS masks each
depth level against the bathymetry, so a level deeper than the shelf turns the
shelf into land.**

One year of weekly releases, 5,300 particles, at four levels:

| level | median lifespan | completed 3 yr | Labrador arrivals | Scotian arrivals | Scotian box masked |
|---|---|---|---|---|---|
| 0.494 m | 155 d | 24% | 52.1% | **1.1%** | 19.4% |
| **47.374 m** | 600 d | 39% | 55.8% | **25.3%** | 29.6% |
| 55.764 m | 555 d | 39% | 55.9% | 19.9% | 32.7% |
| 92.326 m | 1115 d | 54% | 53.5% | 5.6% | 46.3% |

At the surface almost nothing reaches the Scotian Shelf, so the index degenerates
into its Labrador term. At 92 m retention is far better and the shelf is *land*:
57% of the Grand Banks and 46% of the Scotian box are dry, and particles stall
against the mask rather than crossing. The Grand Banks stay navigable to about
40 m and close fast after 65 m, so roughly **47 m** is the only window that is
below the Ekman layer and above the bathymetry.

This is worth stating plainly because a run at the wrong level fails silently: it
produces a clean series of near-zero counts, indistinguishable from a real
absence of transport. `track.py` now prints the masked fraction of both arrival
boxes on every run.

## The result

Five years of weekly releases, 1993–1997, 26,500 particles, at 47.374 m.

Correlation against the published series over the 56-month overlap, on
deseasonalized anomalies — both series carry the annual cycle of the Labrador
Current, and a raw correlation can clear any threshold on shared seasonality
alone:

| | r | p | 95% CI |
|---|---|---|---|
| raw | −0.104 | 0.62 | [−0.38, −0.01] |
| **deseasonalized** | **+0.180** | **0.36** | **[−0.13, +0.28]** |
| detrended | −0.334 | 0.06 | [−0.63, −0.15] |
| both | −0.051 | 0.76 | [−0.40, +0.11] |

Significance is against an autocorrelation-adjusted sample size; consecutive
months are not independent draws, and the usual test on *r* overstates
confidence severalfold here.

## It is not the boxes, and it is not the sample size

Two obvious explanations were tested and both fail.

**The arrival boxes.** They stand in for hydrographic sections whose coordinates
the paper does not publish, and were the prime suspect throughout. Where the
boxes sit does not change the trajectories, only the counting, so 30 placements
were swept over the same 26,500 particles — Labrador moved south away from the
release line, Scotian moved along the shelf. Deseasonalized *r* ranges −0.162 to
**+0.243**, median −0.050, with nothing above 0.4. The surface is smooth and
coherent rather than noisy, and the configured box is already near the best
available. No placement rescues it.

**Counting noise.** These runs release 100 particles a week against the paper's
966, so a noise-dominated index was a live possibility. Splitting each month's
particles into independent halves and comparing the index with itself gives a
split-half agreement of 0.704, a Spearman-Brown reliability of **0.827**, and so
a ceiling of **0.91** on any correlation it could achieve. The index is highly
reproducible. More particles would not help.

The domain was also enlarged from a 45°W eastern edge to 30°W, which raised
retention from 600 to 900 days and moved the index by at most 0.006. Not that
either.

## What this means

The recomputation measures something **stable and repeatable** — reliability
0.83 — that is **not** the published index. That is a more specific finding than
"it didn't work", and it points at the one departure from the paper that could
not be tested: this counts particles that ever *enter a region*, whereas Jutras
et al. count particles *crossing a hydrographic section*.

Those are different quantities, not two estimates of one quantity. A section
crossing is a directed flux at an instant; entering a box is undirected and
cumulative, and counts a particle that drifts in, loiters, and leaves the way it
came. The Labrador box shows the problem in miniature: median time-to-arrival is
30 days, and sweeping it south produces a smooth distance-decay curve rather than
a plateau, so it is largely measuring "drifted a bit south of the release line".

Closing that gap needs the section coordinates, which are not in the paper. They
would have to come from the authors.

## Conclusion

`LCR` cannot be reproduced outside its original framework with what is publicly
available. Monthly fields fail on physics; daily fields with a faithful-looking
Lagrangian pipeline produce a reproducible index that does not track the
published one, and the residual difference is a definition the paper does not
publish.

The pipeline is kept, tested and documented, because it is a few hours from an
answer if those coordinates ever become available. Nothing from it should be
attached as a covariate, and the extension past 2014 should be treated as not
available.

Meanwhile `AMOC` from the RAPID array runs from 2004 to the present and describes
the basin-scale circulation the Labrador Current sits within. It is not a
substitute for `LCR`, but it is available for the years `LCR` is not.

## Reference

Jutras M, Dufour CO, Mucci A, Talbot LC (2023). Large-scale control of the
retroflection of the Labrador Current. *Nature Communications* **14**, 2623.
[doi:10.1038/s41467-023-38321-y](https://doi.org/10.1038/s41467-023-38321-y)
