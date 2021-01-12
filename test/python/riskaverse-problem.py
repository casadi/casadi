from sympy import *
from sympy import Q as QQ

init_printing()

nu = 2
nx = 4
ns = 2
ny = 3
n  = nu + ns + ny
m  = 3

A = MatrixSymbol('A', nx, nx)
Q = MatrixSymbol('Q', nx, nx)
B = MatrixSymbol('B', nx, nu)
R = MatrixSymbol('R', nu, nu)
x0 = MatrixSymbol('x₀', nx, 1)

assumptions.assume.global_assumptions.add(QQ.symmetric(Q))
assumptions.assume.global_assumptions.add(QQ.symmetric(R))

u = MatrixSymbol('u', nu, 1)
s = MatrixSymbol('s', ns, 1)
y = MatrixSymbol('y', ny, 1)
decisionvars = BlockMatrix([[u], [s], [y]])

x1 = A @ x0 + B @ u

scal = lambda s: s * Identity(1)

gs =[
    scal(y[0] - y[1] - s[0]),
    x1.T @ Q @ x1 - scal(s[1]),
    x0.T @ Q @ x0 + u.T @ R @ u - scal(y[0] - y[1] - y[2] - s[1])
]
for i, g in enumerate(gs):
    print(i)
    print(g)
    print("u   : ", simplify(diff(g, u)))
    print("s[0]: ", simplify(diff(g, s[0])))
    print("s[1]: ", simplify(diff(g, s[1])))
    print("y[0]: ", simplify(diff(g, y[0])))
    print("y[1]: ", simplify(diff(g, y[1])))
    print("y[2]: ", simplify(diff(g, y[2])))
