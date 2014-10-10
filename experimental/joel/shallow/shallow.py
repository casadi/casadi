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
from casadi import *
from numpy import *
from matplotlib.pylab import plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import time
import copy

# Physical parameters
g = 9.81 # gravity
poolwidth = 0.2
drag0 = 2.0 # => u(0)
depth0 = 0.01 # => u(1)
sprad = 0.03
spheight = 0.01
endtime = 1.0

# Discretization
numboxes = 20
num_eulersteps = 100
num_measurements = 100

# Plotting
plot_progress = False
numboxes_per_plot = 1

# Discretization
ntimesteps = num_eulersteps*num_measurements
dt = endtime/ntimesteps
dx = poolwidth/numboxes
dy = poolwidth/numboxes
x = linspace(0,poolwidth,numboxes)
y = linspace(0,poolwidth,numboxes)
[X,Y] = meshgrid(x,y)

# Initial conditions
u0 = zeros((numboxes+1, numboxes))
v0 = zeros((numboxes  , numboxes+1))
h0 = zeros((numboxes  , numboxes))
spdist = sqrt( (X-0.04)**2 + (Y-0.04)**2)
I = spdist<sprad/3.0
h0[I] = spheight * cos(3.0*pi*spdist[I]/(2.0*sprad))

# Free parameters
drag = ssym("b")
depth = ssym("H")
p = vertcat([drag,depth])

# The state at a measurement
uk = ssym("uk",numboxes+1, numboxes)
vk = ssym("vk",numboxes  , numboxes+1)
hk = ssym("hk",numboxes  , numboxes)

# Mayer term
hmeas = ssym("h_meas",numboxes  , numboxes)
m = SXFunction([hk,hmeas],[sumAll((hk-hmeas)**2)])
m.init()

# Take one step of the integrator
u = SX(uk)
v = SX(vk)
h = SX(hk)
u[1:-1,:] = u[1:-1,:] + dt*(-g/dx * (h[1:,:]-h[:-1,:]) - u[1:-1,:]*drag)
v[:,1:-1] = v[:,1:-1] + dt*(-g/dy * (h[:,1:]-h[:,:-1]) - v[:,1:-1]*drag)
h[:,:] = h[:,:] + (-depth*dt)*(1.0/dx*(u[1:,:]-u[:-1,:]) + 1.0/dy * (v[:,1:]-v[:,:-1]))

# Create an integrator function
f_step = SXFunction([p,uk,vk,hk],[u,v,h])
f_step.init()

# Expand to SX
#f_step = SXFunction(f_step)
#f_step.init()

# Integrate over one interval
uk = msym("uk",numboxes+1, numboxes)
vk = msym("vk",numboxes  , numboxes+1)
hk = msym("hk",numboxes  , numboxes)
p = msym("p",2)
u,v,h = uk, vk, hk
for j in range(num_eulersteps):
  [u,v,h] = f_step.call([p,u,v,h])

# Create an integrator function
f = MXFunction([p,uk,vk,hk],[u,v,h])
f.init()

print "generated discrete dynamics"

# Allocate memory
u = copy.deepcopy(u0)
v = copy.deepcopy(v0)
h = copy.deepcopy(h0)

# Prepare plotting
if plot_progress:
  fig = plt.figure(1)
  ax = fig.add_subplot(111, projection='3d')
  #plt.clf()
  #plt.grid(True)
  plt.ion()
  plt.hold(False)
  plt.draw()
  plt.show()

# Measurement
h_meas = []

# Simulate once to generate "measurements"
for k in range(num_measurements):
  # Visualize the pool
  if plot_progress:
    #plt.ioff()
    #print h[::numboxes_per_plot]
    ax.cla()
    surf = ax.plot_surface(X[::numboxes_per_plot],Y[::numboxes_per_plot],h[::numboxes_per_plot], rstride=1, cstride=1, cmap=cm.jet, 
      linewidth=0, antialiased=False, vmin = -spheight/10, vmax = spheight/10)
    #plt.contour(X[::numboxes_per_plot],Y[::numboxes_per_plot],h[::numboxes_per_plot])
    #ax.axis([0,poolwidth,0,poolwidth])
    ax.set_zlim(-spheight, spheight)
    plt.draw()
    plt.show()
    #time.sleep(0.02)
    
  # Integrate
  f.setInput([drag0,depth0],0)
  f.setInput(u,1)
  f.setInput(v,2)
  f.setInput(h,3)
  f.evaluate()
  u = f.output(0).toArray()
  v = f.output(1).toArray()
  h = f.output(2).toArray()
    
  # Save a copy of h
  h_meas.append(copy.deepcopy(h))

print "measurements generated"

# Parameters in the single shooting problem
x = msym("x",2)
depth = x[0]
drag = x[1]

# Objective function
obj = 0

# Create expressions for objective functions and constraints
u = MX(u0)
v = MX(v0)
h = MX(h0)
for k in range(num_measurements):
  # Evaluate dynamics
  [u,v,h] = f.call([x,u,v,h])
  
  # Evaluate objective function
  [obj_k] = m.call([h,h_meas[k]])
  
  # add to objective
  obj += obj_k

# Formulate the single shooting NLP
ffcn = MXFunction([x],[obj])
#ffcn.setOption("verbose",True)

ffcn.init()

hfcb = ffcn.jacobian()

#ffcn = SXFunction(ffcn)

#t1 = time.time()
#ffcn.evaluate()
#t2 = time.time()
#print "t = %g" % (t2-t1)


#raise Exception('ss')

# Solve with IPOPT
nlp_solver = IpoptSolver(ffcn)

# Set options
#nlp_solver.setOption("generate_hessian",True)
#nlp_solver.setOption("verbose",True)
nlp_solver.setOption("max_iter",10)
#nlp_solver.setOption("mu_init",1e-10)
nlp_solver.setOption("derivative_test","first-order")

# Initialize NLP solver
nlp_solver.init()

# Set initial condition and bounds
nlp_solver.setInput([drag0,depth0],"x0")
#nlp_solver.setInput([2.0,0.01],"x0")
#nlp_solver.setInput([0.5,0.01],"x0")
nlp_solver.setInput([0,0],"lbx")

# Solve single shooting problem
nlp_solver.evaluate()


