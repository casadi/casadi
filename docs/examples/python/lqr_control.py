#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
from pylab import *
from numpy import *
from casadi import *
from casadi.tools import *
import scipy.linalg
import numpy as np

"""
This example how an Linear Quadratic Regulator (LQR) can be designed
and simulated for a linear-time-invariant (LTI) system
"""

# Problem statement
# -----------------------------------

# System dimensions
ns = 3  # number of states
nu = 2  # number of controls
ny = 2  # number of outputs
N  = 21 # number of control intervals
te = 10 # end time [s]

# The system equations:   x' = A.x + B.u
A = DMatrix([[0.4,0.1,-2],[0,-0.3,4],[1,0,0]])
B = DMatrix([[1,1],[0,1],[1,0]])

# Inspect the open-loop system
[D,V] = linalg.eig(A)
print "Open-loop eigenvalues: ", D

tn = linspace(0,te,N)

# Check if the system is controllable
# -----------------------------------
R = []

p = B
for i in range(ns):
  R.append(p)
  p = mul([A,p])
  
R = horzcat(R)

# R must be of full rank
_, s, _ = linalg.svd(R)
eps = 1e-9
rank =  len([x for x in s if abs(x) > eps])
assert(rank==ns)

# Simulation of the open-loop system
# -----------------------------------

y  = SX.sym("y",ns)
u  = SX.sym("u",nu)

x0 = DMatrix([1,0,0])
# no control
u_ = (DMatrix([[ -1, 1 ],[1,-1]]*((N-1)/2))).T

p = SX.sym("p")

tn = np.linspace(0,te,N)
cdae = SXFunction(controldaeIn(x=y,u=u),[mul(A,y)+mul(B,u)])
cdae.init()


sim = ControlSimulator(cdae,tn)
sim.setOption("integrator", "cvodes")
sim.setOption("integrator_options", {"fsens_err_con": True,"reltol":1e-12})
sim.setOption("nf",20)
sim.init()
sim.setInput(x0,"x0")
sim.setInput(u_,"u")
sim.evaluate()

tf = sim.getMinorT()

figure(1)
plot(tf,sim.getOutput().T)
legend(('s1', 's2','s3'))
title('reference simulation, open-loop, zero controls')
out = sim.getOutput()

# Simulation of the open-loop system
#   sensitivity for initial conditions
# -----------------------------------

x0_pert = DMatrix([0,0,1])*1e-4
sim.setInput(x0+x0_pert,"x0")
sim.setInput(u_,"u")
sim.evaluate()

tf = list(sim.getMinorT())

figure(2)
title('Deviation from reference simulation, with perturbed initial condition')
plot(tf,sim.getOutput().T-out.T,linewidth=3)

# Not supported in current revision, cf. #929
# jacsim = sim.jacobian(CONTROLSIMULATOR_X0,0)
# jacsim.init()
# jacsim.setInput(x0,"x0")
# jacsim.setInput(u_,"u")

# jacsim.evaluate()

# dev_est = []
# for i in range(len(tf)):
#   dev_est.append(mul(jacsim.getOutput()[i*ns:(i+1)*ns,:],x0_pert))

# dev_est = horzcat(dev_est).T
# plot(tf,dev_est,'+k')
# legend(('s1 dev', 's2 dev','s3 dev','s1 dev (est.)', 's2 dev (est.)','s3 dev (est.)'),loc='upper left')


# M = jacsim.getOutput()[-ns:,:]
# # In the case of zero input, we could also use the matrix exponential to obtain sensitivity
# Mref = scipy.linalg.expm(A*te)

# e = fabs(M - Mref)
# e = max(e)/max(fabs(M))
# print e
# assert(e<1e-6)


# Simulation of the open-loop system
#   sensitivity  for controls
# -----------------------------------
# What if we perturb the input?

u_perturb = DMatrix(u_)
u_perturb[0,N/5] = 1e-4
sim.setInput(x0,"x0")
sim.setInput(u_+u_perturb,"u")
sim.evaluate()

figure(3)
title('Deviation from reference simulation, with perturbed controls')
plot(tf,sim.getOutput().T-out.T,linewidth=3)


# Not supported in current revision, cf. #929
# jacsim = sim.jacobian(CONTROLSIMULATOR_U,0)
# jacsim.init()
# jacsim.setInput(x0,"x0")
# jacsim.setInput(u_,"u")

# jacsim.evaluate()

# dev_est = []
# for i in range(len(tf)):
#   dev_est.append(mul(jacsim.getOutput()[i*ns:(i+1)*ns,:],vec(u_perturb)))

# dev_est = horzcat(dev_est).T
# plot(tf,dev_est,'+k')
# legend(('s1 dev', 's2 dev','s3 dev','s1 dev (est.)', 's2 dev (est.)','s3 dev (est.)'),loc='upper left')


# Find feedforward controls
# -----------------------------------
#
#  Our goal is to reach a particular end-state xref_e
#  We can find the necessary controls explicitly in continuous time
#  with the controllability Gramian 
#  
#  http://www-control.eng.cam.ac.uk/jmm/3f2/handout4.pdf
#  http://www.ece.rutgers.edu/~gajic/psfiles/chap5.pdf

x0 = vertcat([1,0,0])
xref_e = vertcat([1,0,0])

states = struct_symSX([
           entry("eAt",shape=(ns,ns)),
           entry("Wt",shape=(ns,ns))
         ])
         
eAt = states["eAt"]
Wt  = states["Wt"]

# We will find the control that allow to reach xref_e at t1
t1 = te

# Initial conditions
e = dense(DMatrix.eye(ns))
states_ = states(0)
states_["eAt"] = e
states_["Wt"] = 0

rhs = struct_SX(states)
rhs["eAt"] = mul(A,eAt)
rhs["Wt"]  = mul([eAt,B,B.T,eAt.T])

dae = SXFunction(daeIn(x=states),daeOut(ode=rhs))
dae.init()

integrator = Integrator("cvodes", dae)
integrator.setOption("tf",t1)
integrator.setOption("reltol",1e-12)
integrator.init()
integrator.setInput(states_,"x0")
integrator.evaluate()

out = states(integrator.getOutput())

Wt_  = out["Wt"]
eAt_ = out["eAt"]

# Check that eAt is indeed the matrix exponential
e = max(fabs(eAt_-scipy.linalg.expm(numpy.array(A*te))))/max(fabs(eAt_))
assert(e<1e-7)

# Simulate with feedforward controls
# -----------------------------------

states = struct_symSX([
          entry("y",shape=ns),      # The regular states of the LTI system
          entry("eAt",shape=(ns,ns))  # The matrix exponential exp(A*(t1-t))
         ])

eAt = states["eAt"]
y   = states["y"]

# Initial conditions
states_ = states(0)
states_["y"] = x0
states_["eAt"] = eAt_

u = mul([B.T,eAt.T,inv(Wt_),xref_e-mul(eAt_,x0)])

rhs = struct_SX(states)
rhs["y"]   = mul(A,y)+mul(B,u)
rhs["eAt"] = -mul(A,eAt)

cdae = SXFunction(controldaeIn(x=states),[rhs])
cdae.init()

# Output function
out = SXFunction(controldaeIn(x=states),[states,u])
out.init()

sim = ControlSimulator(cdae,out,tn)
sim.setOption("integrator", "cvodes")
sim.setOption("integrator_options", {"fsens_err_con": True,"reltol":1e-12})
sim.setOption("nf",20)
sim.init()
sim.setInput(states_,"x0")
sim.evaluate()

e = sim.getOutput()[states.i["y"],-1] - xref_e
assert(max(fabs(e))/max(fabs(xref_e))<1e-6)

tf = sim.getMinorT()
# This will be our reference trajectory for tracking

figure(4)
subplot(211)
title("Feedforward control, states")
plot(tf,sim.getOutput(0)[list(states.i["y"]),:].T)
for i,c in enumerate(['b','g','r']):
  plot(t1,xref_e[i],c+'o')
subplot(212)
title("Control action")
plot(tf,sim.getOutput(1).T)

# Design an infinite horizon LQR
# -----------------------------------

# Weights for the infinite horizon LQR control
Q = DMatrix.eye(ns)
R = DMatrix.eye(nu)

# Continuous Riccati equation
P = SX.sym("P",ns,ns)

ric = (Q + mul(A.T,P) + mul(P,A) - mul([P,B,inv(R),B.T,P]))

dae = SXFunction(daeIn(x=vec(P)),daeOut(ode=vec(ric)))
dae.init()

# We solve the ricatti equation by simulating backwards in time until steady state is reached.
integrator = Integrator("cvodes", dae)
integrator.setOption("reltol",1e-16)
integrator.setOption("stop_at_end",False)
integrator.init()
# Start from P = identity matrix
u = dense(DMatrix.eye(ns))
integrator.reset()
integrator.setInput(vec(u),"x0")
integrator.integrate(0)

# Keep integrating until steady state is reached
for i in range(1,40):
  x0 = integrator.getOutput()
  integrator.integrate(1*i)
  xe = integrator.getOutput()
  e = max(fabs(xe-x0))
  print "it. %02d - deviation from steady state: %.2e" % (i, e)
  if e < 1e-11:
    break
  
# Obtain the solution of the ricatti equation
P_ = integrator.output().reshape((ns,ns))
print "P=", P_

[D,V] = linalg.eig(P_)
assert (min(real(D))>0)
print "(positive definite)"


# Check that it does indeed satisfy the ricatti equation
dae.setInput(integrator.getOutput(),"x")
dae.evaluate()
print max(fabs(dae.getOutput()))
assert(max(fabs(dae.getOutput()))<1e-8)

# From P, obtain a feedback matrix K
K = mul([inv(R),B.T,P_])

print "feedback matrix= ", K

# Inspect closed-loop eigenvalues
[D,V] = linalg.eig(A-mul(B,K))
print "Open-loop eigenvalues: ", D

# Check what happens if we integrate the Riccati equation forward in time
dae = SXFunction(daeIn(x = vec(P)),daeOut(ode=vec(-ric)))
dae.init()

integrator = Integrator("cvodes", dae)
integrator.setOption("reltol",1e-16)
integrator.setOption("stop_at_end",False)
integrator.init()
integrator.setInput(vec(P_),"x0")
integrator.input("x0")[0] += 1e-9 # Put a tiny perturbation
integrator.reset()
integrator.integrate(0)

for i in range(1,10):
  x0 = integrator.getOutput()
  integrator.integrate(0.4*i)
  xe = integrator.getOutput()
  e = max(fabs(xe-x0))
  print "Forward riccati simulation %d; error: %.2e" % (i, e)

# We notice divergence. Why?
stabric = SXFunction([P],[jacobian(-ric,P)])
stabric.init()
stabric.setInput(P_)
stabric.evaluate()

S = stabric.getOutput()

[D,V] = linalg.eig(S)

print "Forward riccati eigenvalues = ", D


# Simulation of the closed-loop system:
#  continuous control action, various continuous references
# ---------------------------------------------------------
#

x0 = DMatrix([1,0,0])

y  = SX.sym("y",ns)

C = DMatrix([[1,0,0],[0,1,0]])
D = DMatrix([[0,0],[0,0]])

temp = inv(blockcat([[A,B],[C,D]]))

F = temp[:ns,-ny:]
Nm = temp[ns:,-ny:]

t = SX.sym("t")

figure(6)

for k,yref in enumerate([ vertcat([-1,sqrt(t)]) , vertcat([-1,-0.5]), vertcat([-1,sin(t)])]):
  u = -mul(K,y) + mul(mul(K,F)+Nm,yref)
  rhs = mul(A,y)+mul(B,u)
  cdae = SXFunction(controldaeIn(t=t, x=y),[rhs])

  # Output function
  out = SXFunction(controldaeIn(t=t, x=y),[y,mul(C,y),u,yref])
  out.init()

  sim = ControlSimulator(cdae,out,tn)
  sim.setOption("integrator", "cvodes")
  sim.setOption("integrator_options", {"fsens_err_con": True,"reltol":1e-12})
  sim.setOption("nf",200)
  sim.init()
  sim.setInput(x0,"x0")
  #sim.setInput(yref_,"u")
  sim.evaluate()

  tf = sim.getMinorT()

  subplot(3,3,1+k*3)
  plot(tf,sim.getOutput(0).T)
  subplot(3,3,2+k*3)
  title('ref ' + str(yref))
  for i,c in enumerate(['b','g']):
    plot(tf,sim.getOutput(1)[i,:].T,c,linewidth=2)
    plot(tf,sim.getOutput(3)[i,:].T,c+'-')
  subplot(3,3,3+k*3)
  plot(tf,sim.getOutput(2).T)


# Simulation of the closed-loop system:
#  continuous control action, continuous feedforward reference
# -----------------------------------------------------------
#  To obtain a continous tracking reference,
#  we augment statespace to construct it on the fly

x0 = vertcat([1,0,0])

# Now simulate with open-loop controls
states = struct_symSX([
           entry("y",shape=ns), # The regular states of the LTI system
           entry("yref",shape=ns), # States that constitute a tracking reference for the LTI system
           entry("eAt",shape=(ns,ns)) # The matrix exponential exp(A*(t1-t))
         ])

y     = states["y"]
eAt   = states["eAt"]

# Initial conditions
states_ = states(0)
states_["y"]    = 2*x0
states_["yref"] = x0
states_["eAt"]  = eAt_


param = struct_symSX([entry("K",shape=(nu,ns))])

param_ = param(0)

uref = mul([B.T,eAt.T,inv(Wt_),xref_e-mul(eAt_,x0)])
u    = uref - mul(param["K"],y-states["yref"])

rhs = struct_SX(states)
rhs["y"]      =  mul(A,y)+mul(B,u)
rhs["yref"]   =  mul(A,states["yref"])+mul(B,uref)
rhs["eAt"]    = -mul(A,eAt)

cdae = SXFunction(controldaeIn(x=states, p=param),[rhs])
cdae.init()

# Output function
out = SXFunction(controldaeIn(x=states, p=param),[states,u,uref,states["yref"]])
out.init()

sim = ControlSimulator(cdae,out,tn)
sim.setOption("integrator", "cvodes")
sim.setOption("integrator_options", {"fsens_err_con": True,"reltol":1e-8})
sim.setOption("nf",20)
sim.init()

# Not supported in current revision, cf. #929
# jacsim = sim.jacobian(CONTROLSIMULATOR_X0,0)
# jacsim.init()

figure(7)
  
for k,(caption,K_) in enumerate([("K: zero",DMatrix.zeros((nu,ns))),("K: LQR",K)]):
  param_["K"] = K_

  sim.setInput(states_,"x0")
  sim.setInput(param_,"p")
  sim.evaluate()
  sim.getOutput()

  tf = sim.getMinorT()
  
  subplot(2,2,2*k+1)
  title('states (%s)' % caption)
  for i,c in enumerate(['b','g','r']):
    plot(tf,sim.getOutput()[states.i["yref",i],:].T,c+'--')
    plot(tf,sim.getOutput()[states.i["y",i],:].T,c,linewidth=2)
  subplot(2,2,2*k+2)
  for i,c in enumerate(['b','g']):
    plot(tf,sim.getOutput(1)[i,:].T,c,linewidth=2)
    plot(tf,sim.getOutput(2)[i,:].T,c+'--')
  title('controls (%s)' % caption)

  # Not supported in current revision, cf. #929
  # # Calculate monodromy matrix
  # jacsim.setInput(states_,"x0")
  # jacsim.setInput(param_,"p")
  # jacsim.evaluate()
  # M = jacsim.getOutput()[-states.size:,:][list(states.i["y"]),list(states.i["y"])]
  
  # # Inspect the eigenvalues of M
  # [D,V] = linalg.eig(M)
  # print "Spectral radius of monodromy (%s): " % caption
  
  # print max(abs(D))

print "Spectral radius of exp((A-BK)*te): "

[D,V] = linalg.eig(scipy.linalg.expm(numpy.array((A-mul(B,K))*te)))
print max(abs(D))

# Simulation of the controller:
# discrete reference, continuous control action
# -----------------------------------------------------------

# Get discrete reference from previous simulation
mi = sim.getMajorIndex()
controls_ = sim.getOutput(2)[:,mi[:-1]]
yref_     = sim.getOutput(3)[:,mi[:-1]]

u_ = vertcat([controls_,yref_])

x0 = DMatrix([1,0,0])

controls = struct_symSX([
             entry("uref",shape=nu),
             entry("yref",shape=ns)
           ])

yref  = SX.sym("yref",ns)
y     = SX.sym("y",ns)
dy    = SX.sym("dy",ns)
u     = controls["uref"]-mul(param["K"],y-controls["yref"])
rhs   = mul(A,y)+mul(B,u)

cdae = SXFunction(controldaeIn(x=y, u=controls, p=param),[rhs])

# Output function
out = SXFunction(controldaeIn(x=y, u=controls, p=param),[y,u,controls["uref"],controls["yref"]])
out.init()

sim = ControlSimulator(cdae,out,tn)
sim.setOption("integrator", "cvodes")
sim.setOption("integrator_options", {"fsens_err_con": True,"reltol":1e-8})
sim.setOption("nf",20)
sim.init()
sim.setInput(2*x0,"x0")
sim.setInput(param_,"p")
sim.setInput(u_,"u")
sim.evaluate()

tf = sim.getMinorT()

figure(8)
subplot(2,1,1)
title('states (%s)' % caption)
for i,c in enumerate(['b','g','r']):
  plot(tf,sim.getOutput(3)[i,:].T,c+'--')
  plot(tf,sim.getOutput()[i,:].T,c,linewidth=2)
subplot(2,1,2)
for i,c in enumerate(['b','g']):
  plot(tf,sim.getOutput(1)[i,:].T,c,linewidth=2)
  plot(tf,sim.getOutput(2)[i,:].T,c+'--')
title('controls (%s)' % caption)

# Not supported in current revision, cf. #929  
# jacsim = sim.jacobian(CONTROLSIMULATOR_X0,0)
# jacsim.init()

# # Calculate monodromy matrix
# jacsim.setInput(x0,"x0")
# jacsim.setInput(param_,"p")
# jacsim.setInput(u_,"u")
# jacsim.evaluate()
# M = jacsim.getOutput()[-ns:,:]

# # Inspect the eigenvalues of M
# [D,V] = linalg.eig(M)
# print "Spectral radius of monodromy, discrete reference, continous control"

# print max(abs(D))
  
# Simulation of the controller:
# discrete reference, discrete control action
# -----------------------------------------------------------

y0     = SX.sym("y0",ns)

u     = controls["uref"]-mul(param["K"],y0-controls["yref"])
rhs   = mul(A,y)+mul(B,u)

cdae = SXFunction(controldaeIn(x=y, x_major=y0, u=controls, p=param),[rhs])

# Output function
out = SXFunction(controldaeIn(x=y, x_major=y0, u=controls, p=param),[y,u,controls["uref"],controls["yref"]])
out.init()

sim = ControlSimulator(cdae,out,tn)
sim.setOption("integrator", "cvodes")
sim.setOption("integrator_options", {"fsens_err_con": True,"reltol":1e-8})
sim.setOption("nf",20)
sim.init()
sim.setInput(2*x0,"x0")
sim.setInput(param_,"p")
sim.setInput(u_,"u")
sim.evaluate()

tf = sim.getMinorT()

figure(9)
subplot(2,1,1)
title('states (%s)' % caption)
for i,c in enumerate(['b','g','r']):
  plot(tf,sim.getOutput(3)[i,:].T,c+'--')
  plot(tf,sim.getOutput()[i,:].T,c,linewidth=2)
subplot(2,1,2)
for i,c in enumerate(['b','g']):
  plot(tf,sim.getOutput(1)[i,:].T,c,linewidth=2)
  plot(tf,sim.getOutput(2)[i,:].T,c+'--')
title('controls (%s)' % caption)

# Not supported in current revision, cf. #929  
# jacsim = sim.jacobian(CONTROLSIMULATOR_X0,0)
# jacsim.init()

# # Calculate monodromy matrix
# jacsim.setInput(x0,"x0")
# jacsim.setInput(param_,"p")
# jacsim.setInput(u_,"u")
# jacsim.evaluate()
# M = jacsim.getOutput()[-ns:,:]

# # Inspect the eigenvalues of M
# [D,V] = linalg.eig(M)
# print "Spectral radius of monodromy, discrete reference, discrete control"

# print max(abs(D))

show()
