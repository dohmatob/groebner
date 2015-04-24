# author: DOHMATOB Elvis Dopgima <gmdopp@gmail.com>

import numbers
from sympy.polys import LT, div
from sympy.polys.polyerrors import ComputationFailed


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


def f4_reduction(F, f):
    """Reduce the multivariate poly modulo the set finite set F of polys"""
    r = 0
    q = [0 for _ in F]
    p = f
    step = 0
    print "F4 reduction of %s mod %s" % (f, F)
    print "+++\r\n"
    while p != 0:
        print "%s: p = %s, q1 = %s, q2 = %s, r = %s" % (step, p, q[0], q[1], r)
        divided = False
        for i, fi in enumerate(F):
            Q, R = div(_LT(p), _LT(fi))
            if R == 0:
                q[i] = (q[i] + Q).expand()
                p = (p - Q * fi).expand()
                divided = True
                break
        if not divided:
            r = (r + _LT(p)).expand()
            p = (p - _LT(p)).expand()
        step += 1
    print "_" * 80
    print "%s = %s + %s" % (
        f, " + ".join(["(%s)(%s)" % (qi, fi) for qi, fi in zip(q, F)
                       if not isinstance(qi, numbers.Number)]), r)
    return q, r
