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

"""
Soft constraints through nlpsol's slack API.

  minimize    f(x, p) + f_s(s, p)
  x, s

  subject to  lbg - S_g s_l <= g(x, p) <= ubg + S_g s_u
              lbx - S_x s_l <=    x    <= ubx + S_x s_u
              0             <=    s    <= ubs      with s = [s_l; s_u]

Two extra entries in the 'nlp' dictionary declare the slack variables:

  's'   the slack variables, 2*ns x 1, stacking [s_l; s_u]
  'f_s' their contribution to the objective; may depend on s and p only

and the 'S' option says which rows they relax: a Sparsity of shape
(ng+nx) x ns holding [S_g; S_x]. It is structural -- its entries count as 1 --
and a structurally empty row leaves that constraint / bound hard.

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
lbg, ubg = [-inf, -inf, -inf], [1.0, 0.5, 3.0]

# Softened rows: g[0], g[2] and the upper bound on x[0]. The empty rows of S
# leave g[1] and the bounds on x[1] hard. S is a Sparsity: only the pattern
# matters, its entries count as 1.
# The norms differ only in how many slacks there are and how they are penalized.
S_SEP = sparsify(DM([[1, 0, 0],    # g[0]  soft
                     [0, 0, 0],    # g[1]  hard
                     [0, 1, 0],    # g[2]  soft
                     [0, 0, 1],    # x[0]  soft
                     [0, 0, 0]])).sparsity()
S_SHARED = sparsify(DM([[1], [0], [1], [1], [0]])).sparsity()

NORMS = {
    # one slack per softened row: sum of violations
    "L1": (S_SEP, lambda s: w * sum1(s)),
    # one slack per softened row: sum of squared violations
    "L2": (S_SEP, lambda s: w * dot(s, s)),
    # a single slack shared by all softened rows: largest violation
    "Linf": (S_SHARED, lambda s: w * sum1(s)),
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

for name, (S, penalty) in NORMS.items():
    ns = S.size2()
    s = SX.sym("s", 2 * ns)

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
