#!/usr/bin/env python3
"""Reference sample points for test/src/unit_tests/cases/pmp_sampling.test.cxx

Independent high-precision implementation of the algorithm in
src/pmp/convert/sample_points.cxx, using Python mpmath (pip install mpmath).

Usage:
  sample_points_reference.py NUM_POINTS BASE [POLE ...]
BASE is an mpmath expression, e.g. 'exp(-1)' or '3-2*sqrt(2)'.
Prints one sample point per line.

Example (prefactor exp(-x) / ((x + 1.5e-10)^2 (x + 1))):
  sample_points_reference.py 10 'exp(-1)' -1.5e-10 -1.5e-10 -1
"""
import sys

import mpmath
from mpmath import mp, mpf, sqrt, log, acos, pi, findroot

BINARY_PRECISION = 3000
OUTPUT_DIGITS = 240
SMALL_POLE_THRESHOLD = mpf('1e-10')


def bracketed_root(f, lo, hi):
    # Illinois method inside the bracket, then polish with the secant method
    x = findroot(f, (lo, hi), solver='illinois', verify=False)
    return findroot(f, x)


def find_b(num_points, base, poles):
    def b_equation(b):
        return (sum(1 - sqrt(-p / (b - p)) for p in poles)
                - b * log(base) / 2 - num_points)

    return bracketed_root(b_equation, SMALL_POLE_THRESHOLD,
                          -2 * num_points / log(base))


def integrated_density(base, poles, b, z):
    # The term for p = 0 is exactly 1 for z > 0
    s = sum(pi if p == 0 else
            2 * acos(sqrt(-p * (b - z) / (b * (z - p))))
            - sqrt(-p / (b - p)) * acos(1 - 2 * z / b) for p in poles)
    return (s - log(base) * (sqrt((b - z) * z)
                             + b / 2 * acos(1 - 2 * z / b))) / pi


def sample_points(num_points, base, poles):
    if num_points == 1:
        return [mpf(0)]
    num_small = min(sum(1 for p in poles if abs(p) <= SMALL_POLE_THRESHOLD),
                    num_points)
    b = find_b(num_points, base, poles)
    points = [None] * num_points
    # Bohr-Sommerfeld roots: integrated_density(z) == n + 1/2
    lo = mpf(0)
    for n in range(num_small, num_points):
        points[n] = bracketed_root(
            lambda z: integrated_density(base, poles, b, z) - n - mpf(1) / 2,
            lo, b)
        lo = points[n]
    # Small points are evenly spaced from 0 to the first root (or to b)
    end = b if num_small == num_points else points[num_small]
    for i in range(num_small):
        points[i] = end * i / num_small
    return points


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    mp.prec = BINARY_PRECISION
    num_points = int(sys.argv[1])
    base = eval(sys.argv[2], {'__builtins__': {}}, vars(mpmath))
    poles = [mpf(p) for p in sys.argv[3:]]
    for x in sample_points(num_points, base, poles):
        print(mp.nstr(x, OUTPUT_DIGITS, min_fixed=-mp.inf, max_fixed=mp.inf)
              if x != 0 else '0')


if __name__ == '__main__':
    main()
