# Licensed to the Apache Software Foundation (ASF) under one or more
# contributor license agreements.  See the NOTICE file distributed with
# this work for additional information regarding copyright ownership.
# The ASF licenses this file to You under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with
# the License.  You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,\
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import argparse
import sys
import os
from mpmath import mp, mpf, mpc, zeta, fsum, floor, __version__ as version
import inspect
import mpmath.functions.zeta as _mpzeta

# Work around mpmath issue #1198 (fixed by PR #1199, unreleased as of Nov-2026):
# _hurwitz_em stops the Euler-Maclaurin tail on an absolute tolerance (|t| < 2^-prec),
# so zeta(s, a) values far below 2^-prec (e.g. ~2^-1022) lose relative accuracy
# (up to ~5e-10 at dps=50). Patch the stopping test to be relative to the result,
# as in PR #1199. Skipped if the installed mpmath already has the fix.
def _patch_hurwitz_em():
  src = inspect.getsource(_mpzeta._hurwitz_em)
  old = "            if ctx.mag(t) < tol:\n"
  if old not in src:
    return False
  src = src.replace(old,
    "            if ctx.mag(t) < tol + ctx.mag(lsum + (-1)**d * tailsum):\n")
  ns = dict(_mpzeta.__dict__)
  exec(src, ns)
  _mpzeta._hurwitz_em = ns['_hurwitz_em']
  return True

patched = _patch_hurwitz_em()

# Negative a (integer s only; otherwise the result is complex).
# mpmath zeta(s, a) for a < 0 is unreliable: below a threshold dps (~2|a|/3 for s=2) it
# returns an mpc with a zero imaginary part and a wrong real value (e.g. zeta(2, -66.26738)
# at dps=30 is -0.01498; true 17.784). Here zeta is only evaluated at positive arguments:
#   reflection: zeta(s, a) = C_s(u) - (-1)^s zeta(s, 1-a), u = a - floor(a),
#               C_s(u) = zeta(s, u) + (-1)^s zeta(s, 1-u). Cost independent of |a|.
#   sum:        zeta(s, a) = sum_{k=0}^{m-1} (a+k)^-s + zeta(s, a+m), m = ceil(-a) + 1.
#               a is a double so each a+k is exact. Cost linear in |a|.
# Each value is computed at dps and 2*dps; the row fails if they do not agree to
# --min-digits significant digits (e.g. cancellation beyond the working precision).
# Note: the reflection identity is the method used by the Java implementation and is for
# development only. Release reference data must be from an independent reference (e.g. the
# sum, or other software).

def zeta_negative(s, a, method):
  s = int(s)
  a = mpf(a)
  if method == 'sum':
    m = int(-floor(a)) + 1
    if m > args.max_sum_terms:
      raise ValueError(f'Sum requires {m} terms > --max-sum-terms {args.max_sum_terms}')
    return fsum((a + k) ** -s for k in range(m)) + zeta(s, a + m)
  u = a - floor(a)
  sign = -1 if s & 1 else 1
  return zeta(s, u) + sign * zeta(s, 1 - u) - sign * zeta(s, 1 - a)

def zeta_checked(s, a, method):
  results = []
  for dps in (args.dps, 2 * args.dps):
    with mp.workdps(dps):
      z = zeta_negative(s, a, method)
      if isinstance(z, mpc):
        raise ArithmeticError(f'Complex result for zeta({s}, {a}) at dps={dps}: {z}')
      results.append(z)
  z1, z2 = results
  with mp.workdps(2 * args.dps):
    if z2 == 0 or abs(z1 - z2) > abs(z2) * mpf(10) ** -args.min_digits:
      raise ArithmeticError(f'zeta({s}, {a}) not converged: dps={args.dps} {z1}; dps={2 * args.dps} {z2}')
  # Round to the output precision
  return +z2

parser = argparse.ArgumentParser(description="Program to compute zeta(s, a).")
parser.add_argument("data", type=str, nargs="+", help="s, a values")
parser.add_argument("--dps", type=int, default=50, help="Digits of precision (default: %(default)s)")
# Default range: all finite doubles
parser.add_argument("--max-a", type=float, default=sys.float_info.max, help="Maximum a in [min, max] (default: %(default)s)")
parser.add_argument("--min-a", type=float, default=-sys.float_info.max, help="Minimum a in [min, max] (default: %(default)s)")
parser.add_argument("--negate-a", default=False, action=argparse.BooleanOptionalAction,
    help="Negate the range a in [min, max] (default: %(default)s)")
parser.add_argument("--negative", choices=['reflection', 'sum'], default='reflection',
    help="Method for a < 0 (default: %(default)s). reflection is for development only")
parser.add_argument("--min-digits", type=int, default=30,
    help="a < 0: digits that must agree at dps and 2*dps (default: %(default)s)")
parser.add_argument("--max-sum-terms", type=int, default=1 << 20,
    help="a < 0: maximum terms for the sum method (default: %(default)s)")
args = parser.parse_args()
mp.dps = args.dps
mp.pretty = True

for filename in args.data:
  sn, n = 0, 0
  with open(filename) as f:
    lines = [line.strip() for line in f]
  def parse(line):
    i = line.index(',') # raise error if not found
    return float(line[0:i]), float(line[i+1:])
  def selected(a):
    return (a > args.max_a or a < args.min_a) == args.negate_a
  data = [parse(line) for line in lines if len(line) and line[0] != '#']
  negative = any(a < 0 and selected(a) for _, a in data)
  with open(f'{os.path.splitext(filename)[0]}.csv', 'w') as f2:
    for line in lines:
      if not len(line) or line[0] == '#':
        print(line, file=f2)
        if 'Generated by' in line:
          print(f'# Evaluated using {os.path.basename(sys.argv[0])}', file=f2)
          print(f'# mpmath ({version}) zeta(s, a) [dps={mp.dps}]', file=f2)
          if patched:
            print('# mpmath patched for issue #1198 (relative Euler-Maclaurin tolerance)', file=f2)
          if negative:
            print(f'# a < 0: {args.negative} with zeta at positive arguments; '
                  f'dps={args.dps} and {2 * args.dps} agree to {args.min_digits} digits', file=f2)
        continue
      s, a = parse(line)
      if not selected(a):
        continue
      if a < 0:
        if s != int(s):
          raise ValueError(f'Non-integer s={s} with a={a} < 0 (complex result)')
        if a == int(a):
          raise ValueError(f'Pole at integer a={a}')
        z = zeta_checked(s, a, args.negative)
      else:
        z = zeta(s, a)
        if isinstance(z, mpc):
          raise ArithmeticError(f'Complex result for zeta({s}, {a}): {z}')
      n += 1
      if abs(float(z)) < 2.2250738585072014e-308:
        sn += 1
      print(f'{s}, {a}, {z}', file=f2)
  if sn:
    print(f'{filename} sub-normal results {sn} / {n} ({100 * sn / n:.3}%)')
