# author: DOHMATOB Elvis Dopgima <gmdopp@gmail.com>

import numbers
import itertools
from sympy import Symbol
from sympy.polys import LT, div, lcm
from sympy.polys.polyerrors import ComputationFailed
from nose.tools import assert_equal

# misc
x, y, z = map(Symbol, "xyz")


def _LT(f):
    try:
        return LT(f)
    except ComputationFailed:
        return f


def _div(a, b):
    try:
        return div(a, b)
    except ComputationFailed:
        return 0, a


def spoly(f, g):
    """Computes the S-polynomial of f and g"""
    lt_f = _LT(f)
    lt_g = _LT(g)
    common = lcm(lt_f, lt_g)
    return ((common * g) / lt_g - (common * f) / lt_f).expand()


def f4_reduction(f, F, verbose=1):
    """Reduce the multivariate poly modulo the finite set polys F

    References
    ----------
    [1] http://algant.eu/documents/theses/thieu.pdf
    """
    r = 0
    q = [0 for _ in F]
    p = f
    step = 0
    if verbose:
        print "F4 reduction of %s mod %s" % (f, F)
        print "+++\r\n"
    while p != 0:
        if verbose:
            print "%s: p = %s, q1 = %s, q2 = %s, r = %s" % (
                step, p, q[0], q[1], r)
        divided = False
        for i, fi in enumerate(F):
            Q, R = _div(_LT(p), _LT(fi))
            if R == 0:
                q[i] = (q[i] + Q).expand()
                p = (p - Q * fi).expand()
                divided = True
                break
        if not divided:
            lt = _LT(p)
            r = (r + lt).expand()
            p = (p - lt).expand()
        step += 1
    if verbose:
        print "_" * 80
        print "%s = %s + %s" % (
            f, " + ".join(["(%s)(%s)" % (qi, fi) for qi, fi in zip(q, F)
                       if not isinstance(qi, numbers.Number)]), r)
    return q, r


def buchberger(F):
    """Bruno Buchberger's algorithm

    Computes Groebner basis for the polynomial ideal spanned by F."""
    G = F
    C = list(itertools.product(F, F))
    cnt = 0
    while C:
        for fi, fj in C:
            C.remove((fi, fj))
            h = f4_reduction(spoly(fi, fj), G, verbose=0)[1]
            print "%3i: spol(%s, %s) mod F = %s" % (cnt, fi, fj, h)
            cnt += 1
            if h != 0:
                C += list(itertools.product(G, [h]))
                G += [h]
    return G


def test_spoly():
    assert_equal(spoly(-1, 3 * x * z ** 2 - y), y)
    assert_equal(spoly(x ** 4 * y - z ** 2, 3 * x * z ** 2 - y),
                 3 * z ** 4 - x ** 3 * y ** 2)
    assert_equal(spoly(4 * x ** 2 * z - 7 * y ** 2,
                       x * y * z ** 2 + 3 * x * z ** 4),
                 12 * x ** 2 * z ** 4 + 7 * y ** 3 * z)


def test_f4_reduction():
    assert_equal(f4_reduction(x ** 2 * y + y - 2, [x * y + 1, x + 1]),
                 ([x, -1], y - 1))
    assert_equal(f4_reduction(x ** 3 + x ** 2 * y + x * y ** 2 + y ** 3,
                              [x * y + 1, x + 1]),
                 ([x + y, x ** 2 - x], y ** 3 - y))


def test_buchberger():
    assert_equal(buchberger([x ** 2 - y * z, y ** 2 - z * x,
                             z ** 2 - x * y]),
                 [x ** 2 - y * z, -x * z + y ** 2, -x * y + z ** 2,
                  y ** 3 - z ** 3])


if __name__ == "__main__":
    print buchberger([x ** 2 - y * z, y ** 2 - z * x, z ** 2 - x * y])
