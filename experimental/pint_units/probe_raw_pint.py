"""Experiments 1+2: raw SX/MX and ca.ArrayInterface inside a pint.Quantity."""
import numpy as np
import casadi as ca, pint, warnings
warnings.simplefilter("ignore")
ureg = pint.UnitRegistry()
def t(name, f):
    try:
        r = f(); print(f"OK   {name}: {type(r).__name__} -> {r!r}")
    except Exception as e:
        print(f"FAIL {name}: {type(e).__name__}: {str(e).splitlines()[0][:150]}")
for label, mk in [("SX", lambda: ca.SX.sym('x')), ("MX", lambda: ca.MX.sym('x')), ("SX vec", lambda: ca.SX.sym('x',3)),
                  ("AI(SX)", lambda: ca.ArrayInterface(ca.SX.sym('x'))), ("AI(MX vec)", lambda: ca.ArrayInterface(ca.MX.sym('x',3)))]:
    print(f"\n=== {label} ===")
    x = mk()
    t("Q(x, m)", lambda: ureg.Quantity(x, 'm'))
    t("x * ureg.m", lambda: x * ureg.m)
    t("ureg.m * x", lambda: ureg.m * x)
    try: q = ureg.Quantity(x, 'm')
    except Exception: continue
    v = ureg.Quantity(x, 'm/s'); tt = ureg.Quantity(2.0, 's')
    t("q + q", lambda: q + q)
    t("q * q", lambda: q*q)
    t("q + 3 cm", lambda: q + 3*ureg.cm)
    t("v * t", lambda: v * tt)
    t("q.to(km)", lambda: q.to('km'))
    t("q ** 2", lambda: q**2)
    t("sin(q/q)", lambda: ca.sin(q/q) if 'AI' not in label else __import__('numpy').sin(q/q))
    t("np.sin(q/m)", lambda: __import__('numpy').sin(q/ureg.m))
    t("np.sqrt(q)", lambda: __import__('numpy').sqrt(q))
    t("q + 1 (dim error expected)", lambda: q + 1)
    t("q.magnitude", lambda: q.magnitude)

print("\n\n######## casadi-namespace functions on raw Quantities ########")
for label, sym in [("SX", ca.SX.sym), ("MX", ca.MX.sym)]:
    print(f"\n=== raw {label} in Quantity ===")
    xs = sym('x', 3); x = xs*ureg.m; v = sym('v')*ureg('m/s')
    t("ca.sin(x) (should be unit error)", lambda: ca.sin(x))
    t("ca.sin(x/ureg.m)", lambda: ca.sin(x/ureg.m))
    t("x[0]", lambda: x[0])
    t("ca.vertcat(x, x)", lambda: ca.vertcat(x, x))
    t("ca.sumsqr(x)", lambda: ca.sumsqr(x))
    t("ca.dot(x,x)", lambda: ca.dot(x, x))
    t("np.sum(x)", lambda: np.sum(x))
    t("ca.jacobian(x**2, x)", lambda: ca.jacobian(x**2, x))
    t("ca.Function([x],[x*v])", lambda: ca.Function('f', [x], [x*v]))
    t("ca.Function([.m],[.m])", lambda: ca.Function('f', [xs, v.magnitude], [(x*v).magnitude]))
    t("x <= 5 m", lambda: x <= 5*ureg.m)
    t("x @ A(1/s)", lambda: x.T @ (np.eye(3)/ureg.s))
    t("ca.mtimes(A, x)", lambda: ca.mtimes(np.eye(3), x))
    t("x.to_base_units", lambda: (x*ureg.km).to_base_units())
