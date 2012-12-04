from casadi import *
t=ssym("t")
x=ssym("x")
rx=ssym("rx")
p=ssym("p")
dp=ssym("dp")

z=ssym("z")
rz=ssym("rz")
rp=ssym("rp")
f = SXFunction(daeIn(**{'x': x, 'z': z}),daeOut(**{'alg': x-z, 'ode': z}))
f.init()
g = SXFunction(rdaeIn(**{'x': x, 'z': z, 'rx': rx, 'rz': rz}),rdaeOut(**{'alg': x-rz, 'ode': rz}))
g.init()

integrator = IdasIntegrator(f,g)
integrator.setOption({'calc_ic': True, 'tf': 2.3, 'reltol': 1e-13, "monitor": ["resB"], "verbose": True, 'augmented_options': {'reltol': 1e-09, 'abstol': 1e-09 }, 'calc_icB': True, 'abstol': 1e-13, 't0': 0.2})
integrator.init()

integrator.input(INTEGRATOR_X0).set(7.1)
if not integrator.input(INTEGRATOR_P).empty():
  integrator.input(INTEGRATOR_P).set(2)
if not integrator.input(INTEGRATOR_RX0).empty():
  integrator.input(INTEGRATOR_RX0).set(0.13)
if not integrator.input(INTEGRATOR_RP).empty():
  integrator.input(INTEGRATOR_RP).set(0.127)

if True:
  integrator.evaluate(0,0)         
  integrator.evaluate(1,0)

integrator.fwdSeed(2).set([1])
integrator.evaluate(1,0) # fail
