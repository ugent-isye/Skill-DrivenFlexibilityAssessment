# Skill-driven flexibility assessment - reference implementation
#
# Copyright (c) 2026, FlandersMake@UGent - ISyE - Ghent University, Belgium
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
Reference implementation of Algorithm 1 of "Hypervolume-based scalar metrics for
skill-driven flexibility assessment in manufacturing system (re)design" (IJPR).

This script is a step-for-step image of that algorithm. See README.md for installation, the
worked example's expected output, and the known limitations. Comments beginning
"Implementation note" cover the practical details the pseudocode leaves open.

Requires the Normaliz executable (https://www.normaliz.uni-osnabrueck.de/).
"""

import atexit
import itertools
import math
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import uuid
from fractions import Fraction

import numpy as np
import polytope as pc


# ---------------------------------------------------------------------------
# Section A -- locating the Normaliz executable
# ---------------------------------------------------------------------------

def _find_normaliz_exe():
    """Locate the Normaliz executable: NORMALIZ_EXE, then PATH, then conda-forge's path."""
    env_path = os.environ.get('NORMALIZ_EXE')
    if env_path and pathlib.Path(env_path).exists():
        return pathlib.Path(env_path)

    on_path = shutil.which('normaliz') or shutil.which('normaliz.exe')
    if on_path:
        return pathlib.Path(on_path)

    candidates = [
        pathlib.Path(sys.prefix) / 'Library' / 'bin' / 'normaliz.exe',  # conda-forge, Windows
        pathlib.Path(sys.prefix) / 'bin' / 'normaliz',                  # conda-forge, Linux/macOS
    ]
    for cand in candidates:
        if cand.exists():
            return cand
    return None


NORMALIZ_EXE = _find_normaliz_exe()


def _normaliz_version():
    if NORMALIZ_EXE is None:
        return 'not found'
    out = subprocess.run([str(NORMALIZ_EXE), '--version'], capture_output=True, text=True,
                          timeout=15).stdout
    m = re.search(r'Normaliz\s+(\S+)', out)
    return m.group(1) if m else 'unknown'


# Versions the expected output was verified against, printed alongside whatever is actually
# installed, so a version difference is visible as a possible cause of any discrepancy.
TESTED_VERSIONS = {
    'python':   '3.11.15',
    'numpy':    '2.4.6',
    'polytope': '0.2.5',
    'normaliz': '3.11.0',
}


def print_version_check():
    installed = {
        'python':   sys.version.split()[0],
        'numpy':    np.__version__,
        'polytope': pc.__version__,
        'normaliz': _normaliz_version(),
    }
    print('Software versions (tested / installed):')
    for lib, tested in TESTED_VERSIONS.items():
        got = installed[lib]
        status = 'ok' if got == tested else 'DIFFERS'
        print(f'  {lib:10s} tested={tested:10s} installed={got:10s} [{status}]')
    print()


# ---------------------------------------------------------------------------
# Section B -- property space, and subsets to polytopes
#                                       [Algorithm 1: step 2, and the P_R / P_Q lines]
#
# A subset is a dict with two keys:
#   'bounds'      -- {property_id: {'min': lo, 'max': hi}} for continuous properties,
#                    {property_id: [allowed_value, ...]} for discrete/categorical ones
#   'constraints' -- list of {'coefficients': {property_id: coefficient}, 'bound': b},
#                    one entry per linear cross-property constraint (g1 or g2)
# ---------------------------------------------------------------------------

def continuous_indices(properties):
    """Indices of the continuous properties (measured with the Lebesgue measure)."""
    return [i for i, p in enumerate(properties) if p['type'] == 'continuous']


def discrete_indices(properties):
    """Indices of the discrete or categorical properties (measured by counting)."""
    return [i for i, p in enumerate(properties) if p['type'] in ('discrete', 'categorical')]


def continuous_bounds_for_subset(sub, properties):
    """Extract {property_id: (lo, hi)} for just the continuous properties of one subset."""
    cidx = continuous_indices(properties)
    bounds = sub.get('bounds', {})
    result = {}
    for gi in cidx:
        pid = properties[gi]['id']
        if pid in bounds:
            b = bounds[pid]
            result[pid] = (b['min'], b['max'])
    return result


def subset_polytope(subset_bounds, subset_constraints, properties):
    """
    Algorithm 1, step 2 / Eq. (A1): the half-space representation P = {x : Ax <= b} of one
    subset's continuous part, given that the discrete/categorical values are already fixed.
    This is the map the appendix writes P_c(.), the subscript naming the discrete and
    categorical combination c at which those values have been fixed.

    Each row is one inequality in the "<=" form Normaliz expects:
      - an upper bound x_i <= hi becomes a 1 at position i, with hi in b;
      - a lower bound x_i >= lo is rewritten as -x_i <= -lo, so a -1 at position i, -lo in b;
      - a cross-property constraint is already in "<=" form, so it is copied as one row.

    `properties` fixes the column order of A. Returns None if the property space has no
    continuous properties at all.
    """
    cidx = continuous_indices(properties)
    n = len(cidx)
    if n == 0:
        return None

    prop_ids = [properties[gi]['id'] for gi in cidx]
    rows_A, rows_b = [], []

    for li, gi in enumerate(cidx):
        pid = properties[gi]['id']
        if pid not in subset_bounds:
            continue
        lo, hi = subset_bounds[pid]

        row_hi = np.zeros(n); row_hi[li] = 1.0
        rows_A.append(row_hi); rows_b.append(hi)

        row_lo = np.zeros(n); row_lo[li] = -1.0
        rows_A.append(row_lo); rows_b.append(-lo)

    for con in (subset_constraints or []):
        row = np.zeros(n)
        if any(pid not in prop_ids for pid in con['coefficients']):
            continue
        for pid, coef in con['coefficients'].items():
            li = prop_ids.index(pid)
            row[li] = coef
        rows_A.append(row)
        rows_b.append(float(con['bound']))

    if not rows_A:
        return None

    A = np.array(rows_A, dtype=float)
    b = np.array(rows_b, dtype=float)
    return pc.Polytope(A, b)


def discrete_combinations(subset_list, properties):
    """
    Algorithm 1: the set K of discrete/categorical combinations admitted by at least one
    subset. Only reachable combinations are enumerated; the rest contribute zero.

    IMPORTANT: pass both families together, as
    `discrete_combinations(cr_subsets + qp_subsets, properties)`. Algorithm 1 uses one
    shared enumeration for all three accumulators.
    """
    didx = discrete_indices(properties)
    if not didx:
        return [()]  # no discrete properties at all: a single "empty" combination

    seen = set()
    result = []
    for sub in subset_list:
        bounds = sub.get('bounds', {})
        per_dim = []
        for gi in didx:
            pid = properties[gi]['id']
            vals = bounds.get(pid, [])
            if not vals:
                per_dim = []
                break
            per_dim.append(sorted(vals))
        if not per_dim:
            continue
        for combo in itertools.product(*per_dim):
            if combo not in seen:
                seen.add(combo)
                result.append(combo)
    return result


def subset_active_at(combo, sub, properties):
    """
    True if this subset is active at the combination `combo`, that is, if its discrete and
    categorical bounds admit it. "Active at" is the appendix's own wording.
    """
    didx = discrete_indices(properties)
    bounds = sub.get('bounds', {})
    for li, gi in enumerate(didx):
        pid = properties[gi]['id']
        if pid in bounds and combo[li] not in bounds[pid]:
            return False
    return True


def polytopes_active_at(combo, subset_list, properties):
    """
    Algorithm 1: the lines building P_R and P_Q, which appendix step 3 names the
    resource-polytope collection and the requirement-polytope collection. Selects the
    subsets admitting `combo`, converts each to its half-space representation, and keeps
    the feasible ones.

    Implementation note: returns None, rather than a list, when the property space has no
    continuous properties at all, so the caller contributes counting measure instead of a
    volume. This is the footnote to appendix step 1.
    """
    if not continuous_indices(properties):
        return None

    polys = []
    for sub in subset_list:
        if not subset_active_at(combo, sub, properties):
            continue
        P = subset_polytope(continuous_bounds_for_subset(sub, properties),
                            sub.get('constraints') or [], properties)
        # { P_c(subset) : subset active at combo }  \  {empty}
        # This test is the appendix's "\ {empty}", not an extra step of our own. is_fulldim
        # solves a Chebyshev-ball linear program: it is True exactly when the inequalities
        # admit a region with interior, which is what carries Lebesgue measure. bool(P) would
        # answer the same question with a Monte Carlo volume estimate, so a thin region could
        # be dropped by chance.
        if P is not None and pc.is_fulldim(P):
            polys.append(P)
    return polys


# ---------------------------------------------------------------------------
# Section C -- volume of one polytope, and of a union   [Algorithm 1: vol(.), UnionVolume]
# ---------------------------------------------------------------------------

_NORMALIZ_WORKDIR = None


def _cleanup_normaliz_workdir():
    """
    Remove the scratch directory when the interpreter exits, but only when it is empty.
    Successful calls delete their own files, so anything left behind belongs to a call that
    failed -- and then the input Normaliz was given is the only record of what was asked of
    it, so the directory is kept and its location reported instead.
    """
    if _NORMALIZ_WORKDIR is None:
        return
    try:
        _NORMALIZ_WORKDIR.rmdir()
    except OSError:
        print(f'normaliz scratch files kept for inspection: {_NORMALIZ_WORKDIR}')


def _normaliz_volume(P):
    """
    Exact volume of a single polytope, via Normaliz, which works directly from the
    half-space representation (P.A, P.b) and needs no vertex enumeration.

    Normaliz picks the algorithm itself, based on the shape of the polytope, and every
    choice it can make is exact (Bruns 2023, "Polytope volume in Normaliz").

    Normaliz is driven as an external process: the polytope is written to a small text
    file, Normaliz is run on it, and its output is parsed back. Two numbers are read: the
    affine dimension, and the "lattice normalized volume" as an exact p/q fraction. That
    lattice volume equals the Euclidean volume times dim!, so dividing by dim! as a
    Fraction recovers the exact Euclidean volume. Converting that Fraction to a float at
    the very end is the ONLY rounding step here. It is deliberately not read from
    Normaliz's own printed "Euclidean volume" line, which is already rounded to about 12
    significant digits.
    """
    global _NORMALIZ_WORKDIR
    if NORMALIZ_EXE is None:
        raise RuntimeError(
            "normaliz executable not found. Set the NORMALIZ_EXE environment variable, "
            "put normaliz on PATH, or install it (e.g. 'conda install -c conda-forge normaliz')."
        )
    if _NORMALIZ_WORKDIR is None:
        _NORMALIZ_WORKDIR = pathlib.Path(tempfile.mkdtemp(prefix='normaliz_'))
        atexit.register(_cleanup_normaliz_workdir)

    A = np.asarray(P.A, dtype=float)
    b = np.asarray(P.b, dtype=float)
    n = A.shape[1]
    name = f'poly_{uuid.uuid4().hex}'
    infile = _NORMALIZ_WORKDIR / f'{name}.in'

    # Normaliz's inhom_inequalities row format is [a_1 ... a_n a_0], meaning a.x + a_0 >= 0.
    # P describes {x : Ax <= b}, i.e. -A_row.x + b >= 0, so each row is -A_row with a_0 = b.
    # float(...) before repr(): NumPy >= 2.0 prints 'np.float64(-1.0)' rather than a bare
    # number, which Normaliz's parser cannot read. Casting loses no precision.
    lines = [f'amb_space {n}', f'inhom_inequalities {A.shape[0]}']
    for row, bi in zip(A, b):
        coeffs = ' '.join(repr(float(-x)) for x in row)
        lines.append(f'{coeffs} {float(bi)!r}')
    lines.append('EuclideanVolume')
    infile.write_text('\n'.join(lines) + '\n')

    subprocess.run([str(NORMALIZ_EXE), '-c', str(infile)],
                    capture_output=True, text=True, cwd=str(_NORMALIZ_WORKDIR), timeout=120)
    outfile = _NORMALIZ_WORKDIR / f'{name}.out'
    if not outfile.exists():
        raise RuntimeError(f'normaliz did not produce output for {infile}')
    out_text = outfile.read_text()

    m_dim = re.search(r'affine dimension of the polyhedron\s*=\s*(\d+)', out_text)
    m_lat = re.search(r'volume \(lattice normalized\)\s*=\s*(\d+)(?:/(\d+))?', out_text)
    m_euc = re.search(r'volume \(Euclidean\)\s*=\s*([0-9.eE+-]+)', out_text)
    if not (m_dim and m_lat and m_euc):
        raise RuntimeError(f'could not find volume/dimension in normaliz output:\n{out_text}')

    dim = int(m_dim.group(1))
    num = int(m_lat.group(1))
    den = int(m_lat.group(2)) if m_lat.group(2) else 1

    exact_volume = Fraction(num, den) / math.factorial(dim)
    result = float(exact_volume)

    # Defensive check only, not the source of the returned value: a parsing bug or a wrong
    # assumption about the dim! normalisation then fails loudly rather than silently.
    reported = float(m_euc.group(1))
    if reported != 0 and abs(result - reported) > 1e-6 * abs(reported):
        raise RuntimeError(f'exact volume {result} disagrees with normaliz-reported {reported}')

    # The .in/.out pair has served its purpose; drop it so that a long run does not
    # accumulate thousands of files. A failing call leaves via one of the raises above and so
    # keeps its files, which is deliberate -- see _cleanup_normaliz_workdir.
    infile.unlink(missing_ok=True)
    outfile.unlink(missing_ok=True)

    return result


def vol(P):
    """
    Algorithm 1's vol(.): the volume of a single polytope.

    Implementation note: a polytope of ambient dimension below two is an interval, so its
    length is read off directly rather than spending a Normaliz call. The result is the same
    either way -- polytope.volume samples inside the bounding box, and in one dimension every
    sample lands inside the polytope.
    """
    if P.dim < 2:
        return P.volume
    return _normaliz_volume(P)


def as_polytopes(R):
    """Normalise a pc.Region, a pc.Polytope, or an empty result to a list of polytopes."""
    if R is None:
        return []
    if isinstance(R, pc.Polytope):
        return [R] if pc.is_fulldim(R) else []
    return list(R.list_poly)


def union_volume(polys):
    """
    Algorithm 1, Function UnionVolume: volume of the union of a collection of polytopes,
    which may overlap. Each polytope after the first contributes only the part of itself
    not already covered, via the sequential set-difference partitioning of Eq. (A3).

    `polys` is the appendix's generic collection P = {H^(1),...,H^(K)}. Step 3 states the
    calculation in that generic form, independently of where the polytopes came from, so the
    same routine serves P_R, P_Q and P_D.

    `as_polytopes(pc.mldivide(P, accumulated))` is the collection the appendix calls Z, the
    disjoint polytopes of H^(k) \ U. The appendix states it as two lines -- form Z, then add
    the volumes of its members -- which is the one expression below.

    The order of the two statements in the loop is load-bearing: the contribution is
    measured against the union BEFORE it is extended. Reversing them would make every
    contribution after the first zero.
    """
    if not polys:
        return 0.0
    v = vol(polys[0])
    accumulated = polys[0]
    for P in polys[1:]:
        v += sum(vol(Z) for Z in as_polytopes(pc.mldivide(P, accumulated)))
        accumulated = pc.Region(as_polytopes(accumulated) + [P])
    return v


def difference_polytopes(polys_Q, polys_R):
    """
    Algorithm 1: the line building P_D, the deficiency-polytope collection of appendix
    step 3, which is the uncovered part of the requirement at one combination.

    The two sides are unioned BEFORE subtracting. A requirement point may be covered
    jointly by two resource subsets and by neither alone, so subtracting subset by subset
    would wrongly report it as uncovered.

    pc.mldivide returns pieces that are already mutually disjoint, so passing them to
    union_volume is correct but does no partitioning work. Appendix step 3 says its generic
    notation H^(k) also covers collections of disjoint polytopes obtained from set
    differences, so P_D goes through UnionVolume like the other two.
    """
    if not polys_Q:
        return []
    if not polys_R:
        # U_Q \ U_R with U_R empty, not a separate case.
        return list(polys_Q)
    return as_polytopes(pc.mldivide(pc.Region(list(polys_Q)), pc.Region(list(polys_R))))


# ---------------------------------------------------------------------------
# Section D -- the metrics                                 [Algorithm 1: main body]
# ---------------------------------------------------------------------------

def _counting_only_totals(cr_subsets, qp_subsets, properties, combinations):
    """
    Implementation note: property spaces with no continuous properties, where each admitted
    combination contributes one unit of counting measure rather than a volume. This is the
    footnote to appendix step 1.
    """
    mu_R = mu_Q = mu_D = 0.0
    for combo in combinations:
        in_R = any(subset_active_at(combo, s, properties) for s in cr_subsets)
        in_Q = any(subset_active_at(combo, s, properties) for s in qp_subsets)
        mu_R += 1.0 if in_R else 0.0
        mu_Q += 1.0 if in_Q else 0.0
        mu_D += 1.0 if (in_Q and not in_R) else 0.0
    return mu_R, mu_Q, mu_D


def compute_metrics(cr_subsets, qp_subsets, properties):
    """
    Algorithm 1 in full: Lambda, Gamma, and Delta_E for one resource against one
    requirement, for one skill.

    Returns the metrics together with the four hypervolumes, under the algorithm's own
    names: mu_R = mu(C_R), mu_Q = mu(Q_P), mu_D = mu(Q_P \\ C_R), mu_I = mu(C_R ∩ Q_P).

    Gamma and Delta_E are nan when mu_R = 0, i.e. when the resource carries no
    functionality for this skill. That is Algorithm 1's own "if mu_R > 0 ... else N/A"
    branch, with nan carrying the N/A; the paper reports those cases as N/A rather than zero
    (appendix step 4). Lambda is 0 when mu_Q = 0, there being nothing to cover.
    """
    combinations = discrete_combinations(cr_subsets + qp_subsets, properties)

    if not continuous_indices(properties):
        mu_R, mu_Q, mu_D = _counting_only_totals(cr_subsets, qp_subsets, properties,
                                                 combinations)
    else:
        mu_R = mu_Q = mu_D = 0.0
        for combo in combinations:
            polys_R = polytopes_active_at(combo, cr_subsets, properties)
            polys_Q = polytopes_active_at(combo, qp_subsets, properties)
            polys_D = difference_polytopes(polys_Q, polys_R)

            mu_R += union_volume(polys_R)
            mu_Q += union_volume(polys_Q)
            mu_D += union_volume(polys_D)

    # Implementation note: guards against a spurious tiny negative from floating-point
    # aggregation when the true intersection is exactly zero.
    mu_I = max(0.0, mu_Q - mu_D)

    # Gamma and Delta_E follow Algorithm 1's "if mu_R > 0 ... else N/A" branch. A resource
    # with no functionality for this skill gives mu_R = 0, and the paper reports those cases
    # as N/A rather than zero (appendix step 4; EL_400 and EL_650 under S2, for example),
    # which is what nan carries here. Algorithm 1 divides mu_I by mu_Q unconditionally; the
    # mu_Q guard below returns 0 when there is nothing to cover.
    Lambda = mu_I / mu_Q if mu_Q > 0 else 0.0
    Gamma = mu_I / mu_R if mu_R > 0 else float('nan')
    Delta_E = mu_D / mu_R if mu_R > 0 else float('nan')

    return dict(mu_R=mu_R, mu_Q=mu_Q, mu_I=mu_I, mu_D=mu_D,
                Lambda=Lambda, Gamma=Gamma, Delta_E=Delta_E)


# ---------------------------------------------------------------------------
# Worked example: NS_600 vs. Scenario 4, skill S1 (cylindrical grinding).
# All values copied directly from Appendix Tables A1 and A3. See README.md.
# ---------------------------------------------------------------------------

PROPERTIES_S1 = [
    {'id': 'p1_diameter',         'type': 'continuous'},
    {'id': 'p2_axial_length',     'type': 'continuous'},
    {'id': 'p3_roughness',        'type': 'continuous'},
    {'id': 'p4_hardness',         'type': 'continuous'},
    {'id': 'p5_location',         'type': 'categorical'},
    {'id': 'p6_workpiece_length', 'type': 'continuous'},
    {'id': 'p7_mass',             'type': 'continuous'},
    {'id': 'p8_profile_tol',      'type': 'continuous'},
]

NS_600_CR = [
    {
        'bounds': {
            'p1_diameter':         {'min': 3,   'max': 380},
            'p2_axial_length':     {'min': 2,   'max': 580},
            'p3_roughness':        {'min': 0.05,'max': 0.8},
            'p4_hardness':         {'min': 30,  'max': 66},
            'p5_location':         ['external'],
            'p6_workpiece_length': {'min': 10,  'max': 600},
            'p7_mass':             {'min': 0.1, 'max': 100},
            'p8_profile_tol':      {'min': 1.7, 'max': 15},
        },
        # g1(x) = rho_m + tau_0 + alpha*x_Lw - x_tprof <= 0, i.e.
        #         alpha*x_Lw - x_tprof <= -(rho_m + tau_0)
        # Table A1: alpha=0.00170, tau_0=1.30, rho_m=0.40  =>  bound = -1.70
        'constraints': [
            {'coefficients': {'p6_workpiece_length': 0.0017, 'p8_profile_tol': -1}, 'bound': -1.7},
        ],
    },
]

SCENARIO4_QP = [
    {
        # external seats
        'bounds': {
            'p1_diameter':         {'min': 40,  'max': 65},
            'p2_axial_length':     {'min': 20,  'max': 60},
            'p3_roughness':        {'min': 0.2, 'max': 0.5},
            'p4_hardness':         {'min': 55,  'max': 60},
            'p5_location':         ['external'],
            'p6_workpiece_length': {'min': 220, 'max': 280},
            'p7_mass':             {'min': 8,   'max': 14},
            'p8_profile_tol':      {'min': 3.5, 'max': 8},
        },
        'constraints': [],
    },
    {
        # internal bore
        'bounds': {
            'p1_diameter':         {'min': 30,  'max': 50},
            'p2_axial_length':     {'min': 20,  'max': 55},
            'p3_roughness':        {'min': 0.3, 'max': 0.7},
            'p4_hardness':         {'min': 55,  'max': 60},
            'p5_location':         ['internal'],
            'p6_workpiece_length': {'min': 220, 'max': 280},
            'p7_mass':             {'min': 8,   'max': 14},
            'p8_profile_tol':      {'min': 4,   'max': 10},
        },
        'constraints': [],
    },
]


EXPECTED = {
    'Lambda': 0.44554455445544555,
    'Gamma': 5.482441377762043e-07,
    'Delta_E': 6.82259371454832e-07,
}

if __name__ == '__main__':
    print_version_check()

    metrics = compute_metrics(NS_600_CR, SCENARIO4_QP, PROPERTIES_S1)
    print('NS_600 vs. Scenario 4 (skill S1, cylindrical grinding)')
    for key in ('mu_R', 'mu_Q', 'mu_I', 'mu_D', 'Lambda', 'Gamma', 'Delta_E'):
        print(f'  {key:15s} = {metrics[key]}')

    for key, expected in EXPECTED.items():
        assert math.isclose(metrics[key], expected, rel_tol=1e-9), (
            f'{key} = {metrics[key]!r} does not match the expected {expected!r}'
        )
    print('\nValues match the expected values.')
