# Numerical reproduction

The original-Hamiltonian DOP853 solver and a separate implicit Radau solver propagate the same two-level problem in different representations. The archive contains their comparison on 91 parameter pairs, followed by an 18-speed calculation of the gain that reaches a specified ratio of excited populations.

The following layout is required only for comparisons with the accepted grid
and figure generation; `numerics/run_point.py` itself can be run without `data/`.

```
data/audited_grid.csv
data/stage3_summary.json
data/calibration_tolerance.csv
figures/
make_manuscript_figures.py
numerics/core.py
numerics/run_point.py
numerics/requirements.txt
```

Python 3.10 or newer is recommended. From the directory containing `numerics/`:

```bash
python -m pip install -r numerics/requirements.txt
python numerics/run_point.py --delta .1 --epsilon .5 --output numerics/local_check.json
python make_manuscript_figures.py
```

All script paths are relative to the downloaded package; no original workspace
path is needed. The scripts can also be invoked by absolute path from another
working directory. The optional output path is relative to the caller's working
directory. Omitting `--output` prints JSON and writes no result file. Figure
generation overwrites only its named figure and diagnostic outputs.

## Quick point check

The default command above uses the accepted grid's reference settings:
`T=8`, boundary order `48`, contour parameter `c=1.7`, local relative tolerance
`3e-14`, and `80` digits for boundary preparation. The expected scaled
probability is

```
C = P exp(pi/delta) = 0.280000506584604  (accepted grid value)
```

Agreement well within `1e-8` relatively is a useful portability check for this
moderate point. It is not an interval error certificate or a substitute for the
broader checks described in the numerical supplement. `example_run.json` is a
fresh run of this one moderate point with the portable entry point; it does not
represent a fresh grid. Typical execution time is seconds on the preparation
machine and is platform dependent.

Direct two-sided matching can be checked by adding `--both-sides`. To vary one
numerical parameter, use e.g. `--T 6`, `--order 40`, `--c 1`, or `--rtol 5e-14`.
These are convergence diagnostics; individual changes need not improve every
digit monotonically.

## Probability grid and independent reproduction

From the package root, run the command below to compare a new Radau solution
with every archived value of the normalized probability F. The last lines
also report the endpoint Fourier integral, one finite-speed example and the
local coefficient in Appendix B, Eq. (B4). The shipped comparison is in
`data/independent_grid.csv`.

```bash
python numerics/independent_check.py --output data/local_comparison.csv
```

The largest relative difference in the shipped comparison was
4.02897e-10 at delta=0.002 and epsilon=0.001. The supplement gives the
contour and accuracy settings. The command does not rerun the original
182 DOP853 propagations.

## Field settings across sweep speeds

Starting at perfect calibration, the script searches outwards in gain error, locates the first sampled bracket for the target and refines the crossing with Brent's method. It checks each inferred setting with the original spinor solver. The article's slow expansion is shown for fixed targets strictly between zero and one. The fast limit follows separately from the pulse area.

~~~bash
python numerics/calibration_tolerance.py --target .01 --output data/local_field_settings.csv
~~~

The 18 rates extend from 0.001 to 100000. The contour offsets and starting times adapt to speed as specified in Supplemental Material Sec. 4.3. At delta = 2 two extra calls propagate the spinor along the real time axis, with the same asymptotic boundary expansion, at the threshold and the slow-limit setting. The function endpoint_quadrature(delta) returns the curvature divided by the bare probability. At delta = 2 their product is 0.302599, giving an amplitude slope of magnitude 0.55009. The columns real_axis_ratio_at_root and real_axis_ratio_at_slow_setting store the results. At rates above 0.1 the inverse_two_term column is empty. Its expansion is used at the slow end of Fig. 1(b).

## One-integral estimate of the tolerance

From the article directory, run the following command after generating the 18 propagated settings:

~~~bash
python numerics/curvature_tolerance.py
~~~

For each speed in data/calibration_tolerance.csv, the script calculates the scaled curvature χ/P_LZ by quadrature and independently from the endpoint propagation at zero gain error. The estimate of the threshold is the square root of the target ratio divided by the scaled curvature. The output data/curvature_tolerance.csv reports the estimated and propagated thresholds, their signed relative difference, and the population ratio obtained by propagating at the estimated threshold. The relative difference uses the propagated threshold as denominator. The plotting script evaluates the same integral separately at 257 logarithmically spaced speeds for the continuous curve.

## Fixed-gain examples

`fixed_gain_examples.py` evaluates the passage at a gain error of 0.1
and sweep rates 0.01, 0.1, and 2. It writes the bare probability
`P_LZ`, the surviving fraction $P_\eta/P_{\mathrm{LZ}}$ from both
solvers, the excited population `P_excited` from Radau, and the
normalized curvature from quadrature.

~~~bash
python numerics/fixed_gain_examples.py --output fixed_gain_check.csv
~~~

## Optional arbitrary-precision check

```bash
python numerics/run_point.py --method mp --delta .1 --epsilon .5 --T 6 --order 80 --dps 55 --degree 48 --c 1 --output numerics/mp_local_check.json
```

The Taylor propagator solves the rational original TDSE independently of
DOP853, but shares its asymptotic boundary recurrence. This check can take much
longer than the binary64 calculation; the full grid should not be scheduled
in arbitrary precision by default. To reproduce the demanding point used in
the numerical audit, change `--delta` and `--epsilon` to `.001`; the other
parameters in the command are the recorded settings. That reproduction can
take several minutes or longer. It was not rerun in preparing this bundle.

The binary64 method's `--dps` controls **only boundary preparation**. Its ODE
propagation remains binary64. The arbitrary-precision method's `--dps` controls
the propagation and its internal tolerance. Both report scaled `C`, avoiding
underflow in the exponentially small probability itself. A small local
tolerance, last boundary term, or inverse-coordinate drift is not a bound on
the global relative probability error.

## Analytical comparisons and endpoints

`run_point.py` first performs propagation, then computes `B0`, `B1`, relative
errors, and any matching accepted-grid discrepancy with `mp.workdps(80)`.
This corrects the old wrapper's default-precision post-processing issue.
The extra precision prevents avoidable rounding in diagnostic arithmetic;
it does not upgrade the accuracy of binary64 input values.

The archived grids use `0 <= epsilon <= 1`, with `eta=1-epsilon`.
The one-point driver also accepts the negative side of the gain error.
At `epsilon=1`, the exact scaled probability is one. At `epsilon=0`, the exact
probability is zero, so a computed nonzero residue is only a numerical floor.
`F=C/epsilon^2` and relative errors are not evaluated at that zero-probability
endpoint. The continuous normalized limit there is a separate linear-response
integral, not a ratio of two floating-point zeros.

## Provenance and limitations

`core.py` extracts `jost_coefficients`, `sval`, `scattering`, and `scattering_mp`
from the previously audited `audit/verify_prefactor.py`. The original wrapper
`correction/numerics.py` only imported these functions. The mathematical
recurrences and propagation are unchanged; unused imports, old batch jobs,
and output-only theoretical comparison fields were removed. Provenance and
the original source hash are in `source_provenance.json`.

The accepted CSV is copied byte-for-byte from the numerical audit. Its SHA256
is recorded in `figure_data_summary.json`. Figure generation uses these data
and the explicit analytical coefficients, with no fitting and no propagation.
The complete grid contains 91 points and originally used 182 propagations.
The archived grid comparisons are empirical diagnostics, not certified
continuum error bounds. Machine and package versions for the new one-point
portability check are in `environment.json`.
