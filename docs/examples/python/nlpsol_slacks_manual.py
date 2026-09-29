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
Soft constraints by hand: the augmentation nlpsol's slack API performs for you.

We solve

  minimize    f(x) + f_s(s)
  x, s

  subject to  lbg - S_g s_l <= g(x) <= ubg + S_g s_u
              lbx - S_x s_l <=  x   <= ubx + S_x s_u
              0             <=  s   <= ubs      with s = [s_l; s_u]

by rewriting it as a plain NLP in the stacked variable z = [x; s_l; s_u]:
every two-sided relaxed row becomes two one-sided rows, and the simple bounds
on x that carry a slack move into the constraint vector.

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
lbg, ubg = [-inf, -inf, -inf], [1.0, 0.5, 3.0]

nx, ng = x.numel(), g.numel()

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

    Sd = DM.ones(S)  # S is structural; its entries count as 1
    S_g = Sd[:ng, :]
    S_x = Sd[ng:, :]

    s = SX.sym("s", 2 * ns)
    s_l, s_u = s[:ns], s[ns:]

    # z = [x; s_l; s_u]
    z = vertcat(x, s)

    # Each relaxed two-sided row splits into two one-sided rows
    G = vertcat(
        g + mtimes(S_g, s_l),  # >= lbg
        g - mtimes(S_g, s_u),  # <= ubg
        x + mtimes(S_x, s_l),  # >= lbx
        x - mtimes(S_x, s_u),  # <= ubx
    )
    lbG = vertcat(DM(lbg), -inf * DM.ones(ng), DM(lbx), -inf * DM.ones(nx))
    ubG = vertcat(inf * DM.ones(ng), DM(ubg), inf * DM.ones(nx), DM(ubx))

    # x is now bounded through G only; the slacks keep their own simple bounds
    lbZ = vertcat(-inf * DM.ones(nx), DM.zeros(2 * ns))
    ubZ = vertcat(inf * DM.ones(nx), inf * DM.ones(2 * ns))

    solver = nlpsol("solver", "sqpmethod", {"x": z, "f": f + penalty(s), "g": G}, opts)
    r = solver(x0=vertcat(DM(x0), DM.zeros(2 * ns)), lbx=lbZ, ubx=ubZ, lbg=lbG, ubg=ubG)

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
