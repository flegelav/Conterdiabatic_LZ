"""Compute the amplitude error at a fixed fraction of the bare excitation."""
import argparse
import csv
import math
from pathlib import Path

from scipy.optimize import brentq

from core import scattering
from independent_check import inverse_settings, scaled_probability


RATES = (.001, .005, .01, .02, .05, .1, .2, .5, 1., 2., 5., 10.,
         30., 100., 300., 1000., 10000., 100000.)


def population_ratio(delta, epsilon):
    if epsilon == 0 or epsilon == 1:
        return float(epsilon)
    return epsilon**2 * scaled_probability(delta, epsilon)


def field_setting(delta, target, leading):
    # Search outwards from perfect calibration to select the inner crossing.
    step = leading / 6
    lo = 0.
    while lo < 1.:
        hi = min(lo + step, 1.)
        if population_ratio(delta, hi) >= target:
            return brentq(lambda eps: population_ratio(delta, eps) - target,
                          lo, hi, xtol=1e-15, rtol=2e-12)
        lo = hi
    raise ArithmeticError("Target could not be bracketed")


def check_ratio(delta, epsilon, real_axis=False):
    # Keep the spinor boundary well into the oscillatory tail as delta grows.
    return float(scattering(delta, epsilon, T=8 * math.sqrt(max(1., delta)),
                            order=48, rtol=3e-14, dps=80,
                            c=min(1.7, math.pi / (4 * delta)),
                            real_axis=real_axis)["C"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=float, default=.01)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    leading, shift = inverse_settings(args.target)
    fields = ["delta", "target_ratio", "leading_epsilon", "inverse_two_term",
              "propagated_epsilon", "dop853_ratio_at_root", "ratio_at_slow_setting",
              "log10_bare_probability", "real_axis_ratio_at_root",
              "real_axis_ratio_at_slow_setting"]
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for delta in RATES:
            root = field_setting(delta, args.target, leading)
            row = dict(delta=delta, target_ratio=args.target, leading_epsilon=leading,
                       propagated_epsilon=root, dop853_ratio_at_root=check_ratio(delta, root),
                       ratio_at_slow_setting=population_ratio(delta, leading),
                       log10_bare_probability=-math.pi / (delta * math.log(10.)))
            if delta <= .1:
                row["inverse_two_term"] = leading - shift * delta**(2/3)
            if delta == 2:
                row["real_axis_ratio_at_root"] = check_ratio(delta, root, real_axis=True)
                row["real_axis_ratio_at_slow_setting"] = check_ratio(delta, leading, real_axis=True)
            writer.writerow(row)
            stream.flush()
            print(delta, root, row["dop853_ratio_at_root"], flush=True)


if __name__ == "__main__":
    main()
