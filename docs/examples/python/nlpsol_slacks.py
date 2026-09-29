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
Soft constraints through nlpsol's slack API.

  minimize    f(x, p) + f_s(s, p)
  x, s

  subject to  lbg - S_lo,g s <= g(x, p) <= ubg + S_up,g s
              lbx - S_lo,x s <=    x    <= ubx + S_up,x s
              0             <=    s    <= ubs

Two extra entries in the 'nlp' dictionary declare the slack variables:

  's'   the slack variables, ns x 1, one per column of S
  'f_s' their contribution to the objective; may depend on s and p only

and the 'S' option says which rows they relax, and on which side: a Sparsity
of shape 2*(ng+nx) x ns stacking [S_lo; S_up], each block [S_g; S_x]. It is
structural -- its entries count as 1 -- and a structurally empty row leaves
that side of that constraint / bound hard.

The solver Function then gains the inputs 's0', 'ubs', 'lam_s0' and the
outputs 's', 'lam_s'; everything else keeps its usual meaning, with 'g' and
'lam_g' reported for the constraints as the user wrote them.

nlpsol_slacks_manual.py performs the same augmentation by hand and prints the
same table.

Joris Gillis, 2026
"""

# ----------------------------------------------------------------------------
# Problem data (identical in nlpsol_slacks_manual.py)
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

# Softened: g[0], g[2] and the simple bound on x[0]; g[1] and x[1] stay hard.
# Rows of [g; x] are numbered 0..n-1; the lower side of row r is row r of S,
# its upper side row n + r.
soft = [0, 2, 3]


def incidence(kind):
    """S for one way of laying the slacks out over the softened rows.

    sym     one slack per row, relaxing both of its sides
    pair    two slacks per row, one per side (independent weights possible)
    up      one slack per row, upper side only; the lower side stays hard
    linf    ONE slack for every softened row and both sides: a symmetric
            worst-case budget, the epigraph form of a max-norm penalty
    """
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
    # sum of violations, either side
    "L1":    (incidence("sym"),  lambda s: w * sum1(s)),
    # sum of squared violations
    "L2":    (incidence("sym"),  lambda s: w * dot(s, s)),
    # lower sides cheap (0.45), upper sides expensive (2.5): both end up used
    "asym":  (incidence("pair"), lambda s: 0.45 * sum1(s[0::2]) + 2.5 * sum1(s[1::2])),
    # only the upper sides may give; pushing a lower bound down is not allowed
    "upper": (incidence("up"),   lambda s: w * sum1(s)),
    # the largest violation over all rows and both sides
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
    s = SX.sym("s", ns)

    nlp = {"x": x, "f": f, "g": g, "s": s, "f_s": penalty(s)}
    solver = nlpsol("solver", "sqpmethod", nlp, dict(opts, S=S))

    r = solver(x0=x0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)

    print("=== %s ===" % name)
    row("f", r["f"])
    row("x", r["x"])
    row("s", r["s"])
    row("g", r["g"])
    row("lam_x", r["lam_x"])
    row("lam_g", r["lam_g"])
    row("lam_s", r["lam_s"])
