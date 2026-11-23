"""Compare curvature estimates of gain tolerance with propagated thresholds."""
import argparse
import csv
import math
from pathlib import Path

from independent_check import endpoint_quadrature, scaled_probability


DATA = Path(__file__).resolve().parent.parent / "data"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--settings", type=Path, default=DATA / "calibration_tolerance.csv")
    parser.add_argument("--output", type=Path, default=DATA / "curvature_tolerance.csv")
    args = parser.parse_args()

    with args.settings.open() as source:
        settings = list(csv.DictReader(source))

    results = []
    for row in settings:
        delta = float(row["delta"])
        target = float(row["target_ratio"])
        exact_setting = float(row["propagated_epsilon"])
        # Both computations return chi/P_LZ directly and work even when P_LZ underflows.
        curvature = endpoint_quadrature(delta)
        curvature_ode = scaled_probability(delta, 0.)
        estimated_setting = math.sqrt(target / curvature)
        achieved_ratio = estimated_setting**2 * scaled_probability(delta, estimated_setting)
        results.append(dict(delta=delta, target_ratio=target, scaled_curvature=curvature,
                            scaled_curvature_ode=curvature_ode,
                            estimated_epsilon=estimated_setting, propagated_epsilon=exact_setting,
                            relative_setting_error=(estimated_setting/exact_setting - 1),
                            achieved_population_ratio=achieved_ratio))

    with args.output.open("w", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    differences = [row["relative_setting_error"] for row in results]
    ode_error = max(abs(row["scaled_curvature_ode"]/row["scaled_curvature"]-1) for row in results)
    print(f"Signed threshold errors over {len(results)} rates: "
          f"{min(differences):.8%} to {max(differences):.8%}")
    print(f"Maximum relative discrepancy in the scaled curvature: {ode_error:.4g}")


if __name__ == "__main__":
    main()
