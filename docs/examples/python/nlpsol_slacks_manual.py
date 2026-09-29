#
#     MIT No Attribution
#
#     Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
from casadi import *
from casadi import *

"""
Soft constraints by hand: the augmentation nlpsol's slack API performs for you.

We solve

  minimize    f(x) + f_s(s)
  x, s

  subject to  lbg - S_lo,g s <= g(x) <= ubg + S_up,g s
              lbx - S_lo,x s <=  x   <= ubx + S_up,x s
              0             <=  s   <= ubs

by rewriting it as a plain NLP in the stacked variable z = [x; s]: every
relaxed two-sided row becomes two one-sided rows, one carrying the lower-side
slacks and one the upper-side slacks, and the simple bounds on x that carry a
slack move into the constraint vector.

Companion of nlpsol_slacks.py, which states the very same problem through the
'S' / 's' / 'f_s' entries of the nlp dictionary. Both print the same table.

Joris Gillis, 2026
"""

# ----------------------------------------------------------------------------
# Problem data (identical in nlpsol_slacks.py)
# ----------------------------------------------------------------------------
w = 1.0

x = SX.sym("x", 2)
f = (x[0] - 3) ** 2 + (x[1] - 2) ** 2
g = vertcat(x[0] + x[1], x[0] - x[1], x[0] + 2 * x[1])

x0 = [0, 0]
lbx, ubx = [-10, -10], [0.5, 10]
lbg, ubg = [-inf, -inf, 4.0], [1.0, 0.5, inf]

nx, ng = x.numel(), g.numel()
n = ng + nx

soft = [0, 2, 3]


def incidence(kind):
    R, C = [], []
    for k, r in enumerate(soft):
        if kind == "sym":
            R += [r, n + r]; C += [k, k]
        elif kind == "pair":
            R += [r, n + r]; C += [2 * k, 2 * k + 1]
        elif kind == "up":
            R += [n + r]; C += [k]
        elif kind == "linf":
            R += [r, n + r]; C += [0, 0]
    return Sparsity.triplet(2 * n, max(C) + 1, R, C)


CASES = {
    "L1":    (incidence("sym"),  lambda s: w * sum1(s)),
    "L2":    (incidence("sym"),  lambda s: w * dot(s, s)),
    "asym":  (incidence("pair"), lambda s: 0.45 * sum1(s[0::2]) + 2.5 * sum1(s[1::2])),
    "upper": (incidence("up"),   lambda s: w * sum1(s)),
    "Linf":  (incidence("linf"), lambda s: w * sum1(s)),
}


def row(label, v):
    v = DM(v)
    print("%-6s%s" % (label, "".join("%14.6f" % float(v[i]) for i in range(v.numel()))))


opts = {
    "print_time": False,
    "print_iteration": False,
    "print_header": False,
    "print_status": False,
    "qpsol": "qrqp",
    "qpsol_options": {"print_iter": False, "print_header": False, "error_on_fail": False},
}

for name, (S, penalty) in CASES.items():
    ns = S.size2()

    Sd = DM.ones(S)  # S is structural; its entries count as 1
    S_lo, S_up = Sd[:n, :], Sd[n:, :]

    s = SX.sym("s", ns)
    z = vertcat(x, s)

    # Each relaxed two-sided row splits into two one-sided rows
    G = vertcat(
        g + mtimes(S_lo[:ng, :], s),  # >= lbg
        g - mtimes(S_up[:ng, :], s),  # <= ubg
        x + mtimes(S_lo[ng:, :], s),  # >= lbx
        x - mtimes(S_up[ng:, :], s),  # <= ubx
    )
    lbG = vertcat(DM(lbg), -inf * DM.ones(ng), DM(lbx), -inf * DM.ones(nx))
    ubG = vertcat(inf * DM.ones(ng), DM(ubg), inf * DM.ones(nx), DM(ubx))

    # x is now bounded through G only; the slacks keep their own simple bounds
    lbZ = vertcat(-inf * DM.ones(nx), DM.zeros(ns))
    ubZ = vertcat(inf * DM.ones(nx), inf * DM.ones(ns))

    solver = nlpsol("solver", "sqpmethod", {"x": z, "f": f + penalty(s), "g": G}, opts)
    r = solver(x0=vertcat(DM(x0), DM.zeros(ns)), lbx=lbZ, ubx=ubZ, lbg=lbG, ubg=ubG)

    # Map the canonical solution back onto the user-facing quantities
    zopt, lamZ, lamG = r["x"], r["lam_x"], r["lam_g"]
    xopt = zopt[:nx]
    sopt = zopt[nx:]
    gopt = Function("g", [x], [g])(xopt)
    lam_s = lamZ[nx:]
    # At most one of the two one-sided rows of a relaxed row can be active
    lam_g = lamG[:ng] + lamG[ng : 2 * ng]
    lam_x = lamZ[:nx] + lamG[2 * ng : 2 * ng + nx] + lamG[2 * ng + nx :]

    print("=== %s ===" % name)
    row("f", r["f"])
    row("x", xopt)
    row("s", sopt)
    row("g", gopt)
    row("lam_x", lam_x)
    row("lam_g", lam_g)
    row("lam_s", lam_s)
