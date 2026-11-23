"""Run one original-Hamiltonian point; numerical propagation has no fitted input.

Example (from any working directory):
  python /path/to/numerics/run_point.py --delta .1 --epsilon .5

The result C=P*exp(pi/delta) stays representable even when P underflows.
Only post-processing uses the two analytical comparison laws.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import time

import mpmath as mp
from core import scattering, scattering_mp

HERE = Path(__file__).resolve().parent


def analysis(result, delta, epsilon):
    """80-digit arithmetic for diagnostics; it does not upgrade numerical data."""
    with mp.workdps(80):
        d, e = mp.mpf(delta), abs(mp.mpf(epsilon))
        C = mp.mpf(str(result["C"]))
        K = mp.pi**2 * mp.mpf(6)**(mp.mpf(2)/3) * mp.gamma(mp.mpf(1)/3) \
            / (10*mp.gamma(mp.mpf(2)/3)**2)
        if e:
            A0 = 2*mp.sin(mp.pi*e/6)
            A1 = K*e*mp.rgamma((1-e)/6)*mp.rgamma((1+e)/6)
            B0, B1 = A0**2/e**2, 2*A0*A1/e**2
            F = C/e**2
        else:
            B0 = mp.pi**2/9
            B1 = 2*mp.pi*K/(3*mp.gamma(mp.mpf(1)/6)**2)
            F = None
        C0, C2 = e**2*B0, e**2*(B0+d**(mp.mpf(2)/3)*B1)
        out = {"arithmetic_dps": 80,
               "precision_note": "These diagnostics retain the accuracy of the numerical input; 80-digit arithmetic does not imply 80-digit accurate propagation.",
               "B0": mp.nstr(B0, 30), "B1": mp.nstr(B1, 30),
               "leading_C": mp.nstr(C0, 30), "two_term_C": mp.nstr(C2, 30)}
        if e and C > 0:
            out.update(F=mp.nstr(F,30),
                       relative_leading_error=mp.nstr(abs(C0-C)/C,25),
                       relative_two_term_error=mp.nstr(abs(C2-C)/C,25),
                       residual_over_delta=mp.nstr((F-B0-d**(mp.mpf(2)/3)*B1)/d,25))
        else:
            out["endpoint_note"] = "At epsilon=0 the exact C is zero; the computed residual is a numerical zero floor. Relative errors and F=C/epsilon^2 are undefined."
        grid = HERE.parent / "data" / "audited_grid.csv"
        if grid.is_file():
            with grid.open() as f:
                for row in csv.DictReader(f):
                    if mp.mpf(row["delta"]) == d and mp.mpf(row["eps"]) == e:
                        ref = mp.mpf(row["C"])
                        out["accepted_grid_C"] = row["C"]
                        out["relative_difference_from_grid"] = mp.nstr(abs(C-ref)/abs(ref),25)
                        break
        return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--delta", default=".1", help="positive adiabatic parameter")
    parser.add_argument("--epsilon", default=".5", help="signed gain error 1-eta in [-1,1]")
    parser.add_argument("--method", choices=["binary64", "mp"], default="binary64")
    parser.add_argument("--T", type=float, help="real endpoint label, default 8 (binary64) or 6 (mp) times sqrt(max(1,delta))")
    parser.add_argument("--order", type=int, help="boundary order, default 48 (binary64), 80 (mp)")
    parser.add_argument("--c", type=float, help="contour offset h=pi/2-c*delta, default 1.7 (binary64) or 1 (mp) when the height is positive")
    parser.add_argument("--dps", type=int, help="boundary digits (binary64) or propagation digits (mp); default 80 or 55")
    parser.add_argument("--rtol", type=float, default=3e-14, help="binary64 local relative tolerance")
    parser.add_argument("--degree", type=int, default=48, help="mp Taylor degree")
    parser.add_argument("--both-sides", action="store_true", help="binary64 direct two-sided matching")
    parser.add_argument("--output", type=Path, help="optional JSON output; otherwise print only")
    args = parser.parse_args()
    d, e = float(args.delta), float(args.epsilon)
    if not (math.isfinite(d) and d > 0 and -1 <= e <= 1):
        parser.error("require a finite delta>0 and -1<=epsilon<=1")
    is_mp = args.method == "mp"
    if args.both_sides and is_mp:
        parser.error("--both-sides is implemented only for binary64")
    T = args.T if args.T is not None else (6 if is_mp else 8) * math.sqrt(max(1., d))
    order = args.order if args.order is not None else (80 if is_mp else 48)
    default_c = 1 if is_mp else 1.7
    if default_c*d >= math.pi/2:
        default_c = math.pi / (4*d)
    c = args.c if args.c is not None else default_c
    dps = args.dps if args.dps is not None else (55 if is_mp else 80)
    if T <= 0 or order < 2 or dps < 25 or c <= 0 or c*d >= float(mp.pi/2):
        parser.error("require T>0, order>=2, dps>=25 and 0<c*delta<pi/2")
    if args.rtol <= 0 or args.degree < 8:
        parser.error("require rtol>0 and degree>=8")
    start = time.perf_counter()
    if is_mp:
        result = scattering_mp(args.delta, args.epsilon, T=T, order=order,
                               dps=dps, degree=args.degree, c=c)
    else:
        result = scattering(d, e, T=T, order=order, rtol=args.rtol, c=c,
                            dps=dps, both_sides=args.both_sides)
    report = {"method": args.method, "result": result,
              "comparison": analysis(result, args.delta, args.epsilon),
              "elapsed_seconds": time.perf_counter()-start,
              "scope": "One-point reproducibility calculation; diagnostics are not certified error bounds."}
    text = json.dumps(report, indent=2)+"\n"
    print(text, end="")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)


if __name__ == "__main__":
    main()
