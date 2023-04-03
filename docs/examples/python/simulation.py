#
#     MIT No Attribution
#
#     Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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
from pylab import *
from scipy.linalg import sqrtm

from casadi import *
from casadi.tools import *

# System states
states = struct_symSX(["x","y","dx","dy"])
x,y,dx,dy = states[...]

# System controls
controls = struct_symSX(["u","v"])
u,v = controls[...]

# System parameters
parameters = struct_symSX(["k","c","beta"])
k,c,beta = parameters[...]

# Provide some numerical values
parameters_ = parameters()
parameters_["k"] = 10
parameters_["beta"] = 1
parameters_["c"] = 1

vel = vertcat(dx,dy)
p = vertcat(x,y)
q = vertcat(u,v)

# System dynamics
F = -k*(p-q) - beta*v*sqrt(sumsqr(vel) + c**2)

# System right hand side
rhs = struct_SX(states)
rhs["x"]  = dx
rhs["y"]  = dy
rhs["dx"] = F[0]
rhs["dy"] = F[1]

# f = SX.fun("f", controldaeIn(x=states,p=parameters,u=controls),daeOut(ode=rhs))


# # Simulation output grid
# N = 100
# tgrid = linspace(0,10.0,N)

# # ControlSimulator will output on each node of the timegrid
# opts = {}
# opts["integrator"] = "cvodes"
# opts["integrator_options"] = {"abstol":1e-10,"reltol":1e-10}
# csim = ControlSimulator("csim", f, tgrid, opts)

# x0 = states(0)

# # Create input profile
# controls_ = controls.repeated(csim.getInput("u"))
# controls_[0,"u"] = 1     # Kick the system with u=1 at the start
# controls_[N/2,"v"] = 2   # Kick the system with v=2 at half the simulation time

# # Pure simulation
# csim.setInput(x0,"x0")
# csim.setInput(parameters_,"p")
# csim.setInput(controls_,"u")
# csim.evaluate()

# output = states.repeated(csim.getOutput())

# # Plot all states
# for k in states.keys():
#   plot(tgrid,output[vertcat,:,k])
# xlabel("t")
# legend(tuple(states.keys()))

# print "xf=", output[-1]

# # The remainder of this file deals with methods to calculate the state covariance matrix as it propagates through the system dynamics

# # === Method 1: integrator sensitivity ===
# # PF = d(I)/d(x0) P0 [d(I)/d(x0)]^T

# P0 = states.squared()
# P0[:,:] = 0.01*DM.eye(states.size)
# P0["x","dy"] = P0["dy","x"] = 0.002

# # Not supported in current revision, cf. #929
# # J = csim.jacobian_old(csim.index_in("x0"),csim.index_out("xf"))
# # J.setInput(x0,"x0")
# # J.setInput(parameters_,"p")
# # J.setInput(controls_,"u")
# # J.evaluate()

# # Jk = states.squared_repeated(J.getOutput())
# # F = Jk[-1]

# # PF_method1 = mtimes([F,P0,F.T])

# # print "State cov (method 1) = ", PF_method1

# # === Method 2: Lyapunov equations ===
# #  P' = A.P + P.A^T
# states_aug = struct_symSX([
#   entry("orig",sym=states),
#   entry("P",shapestruct=(states,states))
# ])

# A = jacobian(rhs,states)

# rhs_aug = struct_SX(states_aug)
# rhs_aug["orig"]  = rhs
# rhs_aug["P"]  = mtimes(A,states_aug["P"]) + mtimes(states_aug["P"],A.T)

# f_aug = SX.fun("f_aug", controldaeIn(x=states_aug,p=parameters,u=controls),daeOut(ode=rhs_aug))

# csim_aug = ControlSimulator("csim_aug", f_aug, tgrid, {"integrator":"cvodes"})

# states_aug(csim_aug.getInput("x0"))["orig"] = x0
# states_aug(csim_aug.getInput("x0"))["P"] = P0

# csim_aug.setInput(parameters_,"p")
# csim_aug.setInput(controls_,"u")
# csim_aug.evaluate()

# output = states_aug.repeated(csim_aug.getOutput())

# PF_method2 = output[-1,"P"]

# print "State cov (method 2) = ", PF_method2

# # === Method 3:  Unscented propagation ===
# # Sample and simulate 2n+1 initial points
# n = states.size

# W0 = 0
# x0 = DM(x0)
# W = DM([  W0 ] + [(1.0-W0)/2/n for j in range(2*n)])

# sqm = sqrtm(n/(1.0-W0)*DM(P0)).real
# sample_x = [ x0 ] + [x0+sqm[:,i] for i in range(n)] + [x0-sqm[:,i] for i in range(n)]

# csim.setInput(parameters_,"p")
# csim.setInput(controls_,"u")

# simulated_x = [] # This can be parallelised
# for x0_ in sample_x:
#   csim.setInput(x0_,"x0")
#   csim.evaluate()
#   simulated_x.append(csim.getOutput()[-1,:])

# simulated_x = vertcat(simulated_x).T

# Xf_mean = mtimes(simulated_x,W)

# x_dev = simulated_x-mtimes(Xf_mean,DM.ones(1,2*n+1))

# PF_method3 = mtimes([x_dev,diag(W),x_dev.T])
# print "State cov (method 3) = ", PF_method3

# show()
