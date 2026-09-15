# Measured: parscale and ranging speedup

Everything in the commits on this branch was justified by counting evaluations in the
code and by tests against stand-in objectives. This is the measurement against FRIDA
itself, which is what those arguments were standing in for.

## How it was run

Two full determinations, launched together on the same idle machine, 6 workers each of
14 cores, `numSample` 100 so the sampling phase does not dominate.

| | branch | `expIDPreString` |
|---|---|---|
| baseline | `145a005` plus **only** the run counter | `PSbase` |
| optimised | this branch | `PSopt` |

`expIDPreString` keeps the two apart everywhere: [config.R:333](config.R#L333) derives
`location.frida` and `location.stella` from `name.output`, so each run gets its own model
copy, simulator copy, worker directories and tmpfs as well as its own output folder.

Both runs report the same starting point: 850 parameters to determine, 37 skipped as
externally ranged, log likelihood in the default case −22627.4440474896.

```
export PATH="$PWD/localSlurm:$PATH"
export FRIDA_COUNT_MODEL_RUNS=TRUE
./submit_UncertaintyAnalysisLevante.sh -s PSopt -w 6 -n 100 -k 10 -a local
```

## Stella runs

The metric the whole plan is built on, and the one that is free of contention between the
two concurrent jobs.

| section | baseline | optimised | saving |
|---|---:|---:|---:|
| parscale determination | 26,289 | 7,994 | **3.29×** |
| range finding | 35,735 | 15,927 | **2.24×** |
| **total** | **62,024** | **23,921** | **2.59×** |

## Wall clock

Determination only, from job start to the counter's report.

| | baseline | optimised | saving |
|---|---:|---:|---:|
| determination | 9h 13m | 4h 36m | **2.00×** |

Wall clock gains less than the run count because not all of the time is Stella: there is R
overhead, file I/O and cluster dispatch per task, and those are unchanged. The two jobs
also overlapped for the first 4h41m and the baseline then had the machine to itself, which
if anything flatters the baseline.

The optimised job completed end to end, determination and sampling, in 4h 41m.

## Does it give the same answer

742 parameters in both, same set and same order.

**Parscale status** — 740 of 742 identical:

| | optimised: determined | notDetermined | skippedExternalRange |
|---|---:|---:|---:|
| baseline: determined | 643 | **2** | 0 |
| notDetermined | 0 | 53 | 0 |
| skippedExternalRange | 0 | 0 | 37 |

Two parameters that the baseline could scale, the optimised run could not. That is the
price of bounding the fallback sweep per parameter rather than globally: the sweep
does not run above the order of the parameter's own range. 2 of 645, 0.3%.

**Parscale values**, where both determined one: 640 of 643 **identical**, 99.5%. The three
that differ do so by at most 1.5e-03 relative.

**Borders**, where both runs determined that border:

| | n | identical | median rel diff | p95 | max |
|---|---:|---:|---:|---:|---:|
| Min | 571 | 95 | 3.86e-09 | 7.84e-07 | 1.87e-02 |
| Max | 585 | 86 | 4.43e-09 | 8.69e-07 | 8.26e-06 |

Nine significant figures at the median. **Only 2 of 571 Min borders differ by more than
1e-3 relative.** The worst is

```
Energy investments.historical cost adjustment[1, Solar, 2]
  baseline   0.1112394557   <- exactly the author's Min bound
  optimised  0.1133183919
  value      0.4638855363,  author range [0.1112394557, 0.8165316169]
```

The baseline border had run all the way into the author bound and been clamped there; the
optimised search stopped just short of it. As a fraction of the sampled half-range the
shift is 5.9e-03, and it is in the conservative direction, a slightly narrower range.

Five Min and two Max borders differ in *whether* they were determined at all.

## What this says about the individual claims

- **The per parameter fallback sweep** was estimated at 2.28× narrower sweeps from a
  synthetic spread. Measured on the real parameter set it is **37 orders down to 14 on
  average, 2.64×** — the run prints this at startup. Its cost is the 2 lost parscales
  above.
- **The parscale determination at 3.29×** is the combined effect of the secant rewrite,
  the per parameter fallback sweep, the parscale cache and the zero-delta probe. It is
  the largest single win and it is nearly free of result changes: 640 of 643 parscales
  identical.
- **`rangeRootTol`** is the only change that moves borders materially, and it is the one
  whose original rationale was already found to be wrong (see
  `developmentTools/testRangeRootTolerance.R`). If exactness matters more than the
  remaining speed, set `rangeRootTol <- NA` in config.R and the border search falls back
  to an absolute 1e-16. This measurement cannot say how much of the range finding's
  2.24× is `rangeRootTol` as against the secant rewrite, the reuse of already evaluated
  points and the single longest-first worker pool; separating them needs another pair of
  runs.

## What is not measured here

- The **cache** work — remembered failed parscales, keyed determinations, reusable
  ranged sampleParms — barely shows: both runs were cold, so nothing was reused. Its
  value is on a *re-run*, where `145a005` already measures 5h44m down to 1h04m, and
  where remembering the failed parscales removes the second pass entirely.
- **Border checks** and the full re-optimisation branch report zero runs, because
  `checkBorderErrors` is FALSE and `treatVarsAsIndep` is TRUE in this config. The
  re-optimisation branch is covered by `developmentTools/testFullReoptBranch.R` instead.
- `numSample` was 100, so nothing here says anything about the sampling phase.
