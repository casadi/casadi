"""Experiment: casadi + pint units.  Run: python demo_casadi_pint.py"""
import warnings
import numpy as np
import casadi as ca
import pint

import casadi_pint as cp
from casadi_pint import CQ

warnings.simplefilter("ignore")      # casadi numpy-mode notice
ureg = pint.UnitRegistry()
CQ.set_registry(ureg)


def check(name, f):
    try:
        print("OK   %-34s %r" % (name, f()))
    except Exception as e:
        print("FAIL %-34s %s: %s" % (name, type(e).__name__, str(e).splitlines()[0][:100]))


print("== 1. CQ arithmetic (pint does the unit algebra) ==")
for kind in (ca.SX, ca.MX):
    x = CQ.sym("x", "km", 3, kind=kind)
    v = CQ.sym("v", "m/s", kind=kind)
    check("x + 5 m", lambda: x + 5 * ureg.m)
    check("x / v  -> s", lambda: (x / v).to("s"))
    check("x * v", lambda: x * v)
    check("x + v (must fail)", lambda: x + v)
    check("x[0] ** 2", lambda: x[0] ** 2)
    check("A @ x  (A in 1/h)", lambda: (np.eye(3) / ureg.h) @ x)
    check("2 m * x  (Quantity on left)", lambda: 2 * ureg.m * x)
    check("5 m + x  (Quantity on left)", lambda: 5 * ureg.m + x)
    check("np.sin(x / 1 m)", lambda: np.sin(x / ureg.m))
    check("cp.sin(x) (must fail)", lambda: cp.sin(x))
    check("cp.sqrt(x*v)", lambda: cp.sqrt(x * v))
    check("cp.vertcat(x, 5 m)", lambda: cp.vertcat(x, 5 * ureg.m))
    check("cp.jacobian(x*v, v)", lambda: cp.jacobian(x * v, v))
    check("x <= 90 km/h * 1 min", lambda: x <= 90 * ureg("km/h") * ureg.min)

print("\n== 2. passing CQ into plain casadi via __SX__/__MX__ ==")
x = CQ.sym("x", "m")
check("ca.SX(x)", lambda: ca.SX(x))
check("ca.MX(x) (SX-backed -> refuse)", lambda: ca.MX(x))
check("ca.vertcat(x, x)", lambda: ca.vertcat(x, x))
check("ca.Function('f',[x],[x*x])", lambda: ca.Function("f", [x], [x * x]))
check("ca.sin(x/(1 cm)) -> 100 folded", lambda: ca.sin(x / ureg.cm))
check("ca.sin(x)  (non-strict: silent!)", lambda: ca.sin(x))
CQ.strict = True
check("ca.sin(x)  (strict)", lambda: ca.sin(x))
check("ca.sin(x/(1 cm))  (strict)", lambda: ca.sin(x / ureg.cm))
CQ.strict = False
check("float(CQ(2 km / 4 m))", lambda: float(CQ(2 * ureg.km / (4 * ureg.m))))

print("\n== 3. unit-aware Function ==")
x = CQ.sym("x", "km")
v = CQ.sym("v", "km/h")
f, iu, ou = cp.Function("travel_time", [x, v], [(x / v).to("min")])
print(f, iu, ou)
print("f(3 km, 60 km/h) =", f(3, 60) * ou[0])

print("\n== 4. small optimal-control problem in mixed units (Opti) ==")
# Drive D = 1.2 km in T = 1.5 min, speed <= 72 km/h, |accel| <= 2 m/s^2,
# start and end at rest; minimise integral of accel^2.
N = 30
opti = ca.Opti()
dt = (1.5 * ureg.min / N)
v = CQ(opti.variable(N + 1), "km/h")          # decision vars in km/h
a = CQ(opti.variable(N), "m/s**2")            # ... and m/s^2
for k in range(N):
    opti.subject_to(v[k + 1] == v[k] + a[k] * dt)   # pint converts both sides
pos = cp.sum1((v[0:N] + v[1:N + 1]) / 2) * dt
opti.subject_to(pos == 1.2 * ureg.km)
opti.subject_to(v[0] == 0 * ureg("km/h"))
opti.subject_to(v[N] == 0 * ureg("km/h"))
opti.subject_to(opti.bounded(0, v.m_as("km/h"), 72))
opti.subject_to(cp.fabs(a) <= 2 * ureg("m/s**2"))
J = cp.sumsqr(a) * dt
print("objective unit:", J.units)
opti.minimize(J.m_as("m**2/s**3"))
opti.solver("ipopt", {"print_time": False}, {"print_level": 0, "sb": "yes"})
sol = opti.solve()
vs = sol.value(v.expr)
print("J = %.4f m^2/s^3" % sol.value(J.m_as("m**2/s**3")))
print("v_max = %.2f km/h, distance = %.4f km"
      % (vs.max(), sol.value(pos.m_as("km"))))
