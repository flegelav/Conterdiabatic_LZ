# Calibration tolerance for a counterdiabatic Landau–Zener sweep

This repository contains the numerical data, solvers, and figures for the
study of gain errors in a counterdiabatic Landau–Zener sweep. The gain is
$\eta$, the applied pulse amplitude divided by the ideal amplitude, and
the signed gain error is $\varepsilon=1-\eta$. The sweep rate $\delta$ is
dimensionless. The surviving fraction is $R=P_\eta/P_{\mathrm{LZ}}$, with
$P_{\mathrm{LZ}}=\exp(-\pi/\delta)$.

| File | Contents |
| --- | --- |
| [audited_grid.csv](data/audited_grid.csv) | Transition probabilities and asymptotic comparisons for 91 parameter pairs. `C` is $R$, and `F` is $R/\varepsilon^2$. |
| [independent_grid.csv](data/independent_grid.csv) | Comparison of the archived probabilities with separate Radau integrations on the same grid. |
| [calibration_tolerance.csv](data/calibration_tolerance.csv) | Gain-error thresholds for a hundredfold reduction of excitation at 18 sweep rates, with cross-checks. |
| [curvature_tolerance.csv](data/curvature_tolerance.csv) | Threshold estimates from the curvature, with the actual excitation ratios at the estimated gains. |
| [fixed_gain_examples.csv](data/fixed_gain_examples.csv) | Three calculations at $\varepsilon=0.1$, including probabilities and results from both solvers. |
| [stage3_summary.json](data/stage3_summary.json) | Numerical checks and aggregate errors for the 91-point grid. |
| [figures/](figures/) | Main and supplemental figures in vector PDF and PNG. |

## Reproduce the calculations

Use Python 3.10 or newer. Run the commands from the repository root.

~~~bash
python -m pip install -r numerics/requirements.txt
python numerics/run_point.py --delta 0.1 --epsilon 0.5
python make_manuscript_figures.py
~~~

The point calculation gives $R\approx0.280000506584604$. The plotting
script reads the archived CSV files and updates the three figures and
`figure_data_summary.json`.

To regenerate the tolerance tables, compare the archived grid with the
Radau solver, or reproduce the fixed-gain examples, write to a separate
directory.

~~~bash
mkdir -p reproduced
python numerics/calibration_tolerance.py --target 0.01 --output reproduced/calibration_tolerance.csv
python numerics/curvature_tolerance.py --settings reproduced/calibration_tolerance.csv --output reproduced/curvature_tolerance.csv
python numerics/independent_check.py --output reproduced/independent_grid.csv
python numerics/fixed_gain_examples.py --output reproduced/fixed_gain_examples.csv
~~~

See the [numerical documentation](numerics/README.md) for solver settings,
column definitions, and the arbitrary-precision check. The
[manifest](manifest.json) records SHA-256 digests for the archived files.
