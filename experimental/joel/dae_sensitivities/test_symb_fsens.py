from casadi import *
from copy import deepcopy
#deepcopy = DMatrix

print "Testing sensitivity analysis in CasADi"

# Time 
t = ssym("t")

# Differential states
s = ssym("s"); v = ssym("v"); m = ssym("m")
x = vertcat([s,v,m])

# State derivatives
sdot = ssym("sdot"); vdot = ssym("vdot"); mdot = ssym("mdot")
xdot = vertcat([sdot,vdot,mdot])

# Control
u = ssym("u")

# Constants
alpha = 0.05 # friction
beta = 0.1   # fuel consumption rate
  
# Differential equation
ode = vertcat([
  v-sdot,
  (u-alpha*v*v)/m - vdot,
  -beta*u*u       - mdot])
  
# Quadrature
quad = v**3 + ((3-sin(t)) - u)**2

# DAE callback function
ffcn = SXFunction(daeIn(t=t,x=x,xdot=xdot,p=u),daeOut(ode=ode,quad=quad))

# Time length
tf = 0.5

# Initial position
x0 = [0.,0.,1.]

# Parameter guess
u0 = 0.4

# Integrator
for MyIntegrator in (CVodesIntegrator,CollocationIntegrator):
  print "========"
  print "Integrator: ", MyIntegrator.__name__
  print "========"

  # Integrator options
  opts = {}
  opts["tf"]=tf
  if MyIntegrator==CollocationIntegrator:
    opts["implicit_solver"] = KinsolSolver
    opts["implicit_solver_options"] = {"linear_solver_creator":CSparse}

  # Integrator
  I = MyIntegrator(ffcn)
  I.setOption(opts)
  I.init()

  # Integrate to get results
  I.setInput(x0,INTEGRATOR_X0)
  I.setInput(u0,INTEGRATOR_P)
  I.evaluate()
  xf = deepcopy(I.output(INTEGRATOR_XF))
  qf = deepcopy(I.output(INTEGRATOR_QF))
  print "%50s" % "Unperturbed solution:", "xf  = ", xf, ", qf  = ", qf

  # Perturb solution to get a finite difference approximation
  h = 0.001
  I.setInput(u0+h,INTEGRATOR_P)
  I.evaluate()
  xf_pert = deepcopy(I.output(INTEGRATOR_XF))
  qf_pert = deepcopy(I.output(INTEGRATOR_QF))
  print "%50s" % "Finite difference approximation:", "d(xf)/d(p) = ", (xf_pert-xf)/h, ", d(qf)/d(p) = ", (qf_pert-qf)/h

  # Calculate once, forward
  I_fwd = I.derivative(1,0)
  I_fwd.setInput(x0,INTEGRATOR_X0)
  I_fwd.setInput(u0,INTEGRATOR_P)
  I_fwd.setInput([0,0,0],INTEGRATOR_NUM_IN+INTEGRATOR_X0)
  I_fwd.setInput(1.0,INTEGRATOR_NUM_IN+INTEGRATOR_P)
  I_fwd.evaluate()
  fwd_xf = deepcopy(I_fwd.output(INTEGRATOR_NUM_OUT+INTEGRATOR_XF))
  fwd_qf = deepcopy(I_fwd.output(INTEGRATOR_NUM_OUT+INTEGRATOR_QF))
  print "%50s" % "Forward sensitivities:", "d(xf)/d(p) = ", fwd_xf, ", d(qf)/d(p) = ", fwd_qf

  # Calculate once, adjoint
  I_adj = I.derivative(0,1)
  I_adj.setInput(x0,INTEGRATOR_X0)
  I_adj.setInput(u0,INTEGRATOR_P)
  I_adj.setInput([0,0,0],INTEGRATOR_NUM_IN+INTEGRATOR_XF)
  I_adj.setInput(1.0,INTEGRATOR_NUM_IN+INTEGRATOR_QF)
  I_adj.evaluate()
  adj_x0 = deepcopy(I_adj.output(INTEGRATOR_NUM_OUT+INTEGRATOR_X0))
  adj_p = deepcopy(I_adj.output(INTEGRATOR_NUM_OUT+INTEGRATOR_P))
  print "%50s" % "Adjoint sensitivities:", "d(qf)/d(x0) = ", adj_x0, ", d(qf)/d(p) = ", adj_p

  # Perturb adjoint solution to get a finite difference approximation of the second order sensitivities
  I_adj.setInput(x0,INTEGRATOR_X0)
  I_adj.setInput(u0+h,INTEGRATOR_P)
  I_adj.setInput([0,0,0],INTEGRATOR_NUM_IN+INTEGRATOR_XF)
  I_adj.setInput(1.0,INTEGRATOR_NUM_IN+INTEGRATOR_QF)
  I_adj.evaluate()
  adj_x0_pert = deepcopy(I_adj.output(INTEGRATOR_NUM_OUT+INTEGRATOR_X0))
  adj_p_pert = deepcopy(I_adj.output(INTEGRATOR_NUM_OUT+INTEGRATOR_P))
  print "%50s" % "FD of adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", (adj_x0_pert-adj_x0)/h, ", d2(qf)/d(p)d(p) = ", (adj_p_pert-adj_p)/h

  # Forward over adjoint to get the second order sensitivities
  I_adj.setInput(x0,INTEGRATOR_X0)
  I_adj.setInput(u0,INTEGRATOR_P)
  I_adj.setFwdSeed(1.0,INTEGRATOR_P)
  I_adj.setInput([0,0,0],INTEGRATOR_NUM_IN+INTEGRATOR_XF)
  I_adj.setInput(1.0,INTEGRATOR_NUM_IN+INTEGRATOR_QF)
  I_adj.evaluate(1,0)
  fwd_adj_x0 = deepcopy(I_adj.fwdSens(INTEGRATOR_NUM_OUT+INTEGRATOR_X0))
  fwd_adj_p = deepcopy(I_adj.fwdSens(INTEGRATOR_NUM_OUT+INTEGRATOR_P))
  print "%50s" % "Forward over adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", fwd_adj_x0, ", d2(qf)/d(p)d(p) = ", fwd_adj_p

  # Adjoint over adjoint to get the second order sensitivities
  I_adj.setInput(x0,INTEGRATOR_X0)
  I_adj.setInput(u0,INTEGRATOR_P)
  I_adj.setInput([0,0,0],INTEGRATOR_NUM_IN+INTEGRATOR_XF)
  I_adj.setInput(1.0,INTEGRATOR_NUM_IN+INTEGRATOR_QF)
  I_adj.setAdjSeed(1.0,INTEGRATOR_NUM_OUT+INTEGRATOR_P)
  I_adj.evaluate(0,1)
  adj_adj_x0 = deepcopy(I_adj.adjSens(INTEGRATOR_X0))
  adj_adj_p = deepcopy(I_adj.adjSens(INTEGRATOR_P))
  print "%50s" % "Adjoint over adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", adj_adj_x0, ", d2(qf)/d(p)d(p) = ", adj_adj_p

  # Operator overloading approach
  I.setInput(x0,INTEGRATOR_X0)
  I.setInput(u0,INTEGRATOR_P)
  I.setFwdSeed([0,0,0],INTEGRATOR_X0)
  I.setFwdSeed(1.0,INTEGRATOR_P)
  I.reset(1,0,0)
  I.integrate(tf)
  oo_xf = deepcopy(I.fwdSens(INTEGRATOR_XF))
  oo_qf = deepcopy(I.fwdSens(INTEGRATOR_QF))
  print "%50s" % "Forward sensitivities via OO:", "d(xf)/d(p) = ", oo_xf, ", d(qf)/d(p) = ", oo_qf
  