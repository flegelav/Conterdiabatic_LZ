"""Reproduce the three tabulated examples at a gain error of 0.1."""
import argparse
import csv
import math
from pathlib import Path

from calibration_tolerance import check_ratio
from independent_check import endpoint_quadrature, scaled_probability


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="output CSV path")
    args = parser.parse_args()

    epsilon = 0.1
    rows = []
    for delta in (0.01, 0.1, 2.0):
        bare_probability = math.exp(-math.pi / delta)
        scaled_ratio = scaled_probability(delta, epsilon)
        radau_ratio = epsilon**2 * scaled_ratio
        rows.append({
            "delta": delta,
            "epsilon": epsilon,
            "eta": 1 - epsilon,
            "P_LZ": bare_probability,
            "R_Radau": radau_ratio,
            "R_DOP853": check_ratio(delta, epsilon),
            "P_excited": bare_probability * radau_ratio,
            "R_over_epsilon_squared": scaled_ratio,
            "chi_over_P_LZ": endpoint_quadrature(delta),
        })

    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
