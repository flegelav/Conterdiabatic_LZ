"""Direct channel propagation for a fixed-shape counterdiabatic Landau–Zener sweep."""
import argparse
import csv
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import gamma, rgamma


def coefficients(eps):
    """Leading and first subleading coefficients of the scaled probability."""
    # These are A_0/epsilon and A_1/epsilon, continued to epsilon=0.
    leading_amplitude = np.pi / 3 * np.sinc(eps / 6)
    correction_amplitude = (
        (9/2)**(1/3) * np.pi**2 * gamma(1/3) / (5 * gamma(2/3)**2)
        * rgamma((1 + eps)/6) * rgamma((1 - eps)/6)
    )
    return leading_amplitude**2, 2 * leading_amplitude * correction_amplitude


def inverse_settings(target):
    """Return the leading setting and its positive finite-rate shift coefficient."""
    if not 0 < target < 1:
        raise ValueError("the target fraction must be strictly between zero and one")
    angle = math.asin(math.sqrt(target) / 2)
    leading = 6 * angle / math.pi
    b0, b1 = coefficients(leading)
    shift = 3 * math.tan(angle) * (b1 / b0) / math.pi
    return leading, shift


def P_lz_cd(delta, eta, *, logarithm=False):
    """Two-term, small-rate probability; logarithm=True returns log(P)."""
    if not math.isfinite(delta) or delta <= 0 or not 0 <= eta <= 2:
        raise ValueError("delta must be finite and positive and eta must be in [0, 2]")
    eps = abs(1 - eta)
    if eps == 0:
        return -math.inf if logarithm else 0.0
    b0, b1 = coefficients(eps)
    log_p = -np.pi/delta + 2*math.log(eps) + math.log(b0 + delta**(2/3)*b1)
    return log_p if logarithm else math.exp(log_p)


def contour(delta, offset):
    k = 1 / (2 * delta)
    if not 0 < offset < k*np.pi/2:
        raise ValueError("the contour requires 0 < offset < pi/(4*delta)")
    angle = brentq(lambda a: k*(np.pi/2-a-np.sin(2*a)/2)-offset, 0, np.pi/2)
    return k, angle


def scaled_probability(delta, eps, *, offset=None, cutoff=None, tol=2e-11):
    if offset is None:
        offset = min(.75, np.pi/(8*delta))
    if cutoff is None:
        cutoff = 8 + max(0., .5*math.log(delta))
    k, angle = contour(delta, offset)

    def derivative(t, state):
        a, b = np.ascontiguousarray(state).view(np.complex128)
        h = np.cosh(t + 1j*angle)
        g = 0.5j/h
        return np.array([-4j*k*h*h*a + g*b, eps*eps*g*a]).view(float)

    start = np.array([0j, 1+0j]).view(float)
    sol = solve_ivp(derivative, (-cutoff, 0), start, method='Radau',
                    rtol=tol, atol=tol/100)
    if not sol.success:
        raise RuntimeError(sol.message)
    a, b = np.ascontiguousarray(sol.y[:, -1]).view(np.complex128)
    return (2*math.exp(2*offset)*(a*b.conjugate()).imag)**2


def endpoint_quadrature(delta, offset=None):
    """Return chi divided by the bare Landau–Zener probability."""
    if offset is None:
        offset = min(.75, np.pi/(8*delta))
    cutoff = 8 + max(0., .5*math.log(delta))
    k, angle = contour(delta, offset)

    def integrand(t):
        u = t + 1j*angle
        s = u + 0.5*np.sinh(2*u)
        return (np.exp(2*k*(np.pi/2 + 1j*s))/np.cosh(u)).real

    amplitude, err = quad(integrand, 0, cutoff, epsabs=3e-13, epsrel=3e-13, limit=400)
    return amplitude**2


def local_coefficient(r):
    """Extract the first subleading coefficient directly on the physical axis."""
    with mp.workdps(90):
        r = mp.mpf(r)
        d = mp.betainc(mp.mpf(3)/2, mp.mpf(1)/2, 0, r)
        return d**(mp.mpf(1)/3)*(1/(2*r**(mp.mpf(3)/2))-1/(3*d))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write the full grid comparison as CSV")
    args = parser.parse_args()
    grid = Path(__file__).resolve().parent.parent / "data" / "audited_grid.csv"
    with grid.open() as stream:
        points = list(csv.DictReader(stream))
    rows = []
    for original in points:
        delta, eps = float(original["delta"]), float(original["eps"])
        reference = float(original["F"])
        actual = scaled_probability(delta, eps)
        rows.append(dict(delta=delta, epsilon=eps, F_archive=reference,
                         F_independent=actual, relative_difference=abs(actual/reference-1)))
    if args.output:
        with args.output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0])
            writer.writeheader()
            writer.writerows(rows)
    worst = max(rows, key=lambda r: r["relative_difference"])
    print("Number of parameter pairs:", len(rows))
    print("Largest relative difference in F:", worst)
    print("Endpoint F(0.01, 0) by direct quadrature:", endpoint_quadrature(.01))
    print("F(0.01, 0.1) by direct propagation:", scaled_probability(.01, .1))
    print("Estimated a0 at r=1e-20:", mp.nstr(local_coefficient("1e-20"), 25))


if __name__ == "__main__":
    main()
