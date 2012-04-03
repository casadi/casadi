from casadi import *
from numpy import *
from matplotlib.pylab import plt
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
numboxes = 30
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
hmeas = ssym("h_meas",numboxes  , numboxes)

# Integrate over one interval with Euler forward
u = SXMatrix(uk)
v = SXMatrix(vk)
h = SXMatrix(hk)
for j in range(num_eulersteps):
  u[1:-1,:] = u[1:-1,:] + dt*(-g/dx * (h[1:,:]-h[:-1,:]) - u[1:-1,:]*drag)
  v[:,1:-1] = v[:,1:-1] + dt*(-g/dy * (h[:,1:]-h[:,:-1]) - v[:,1:-1]*drag)
  h[:,:] = h[:,:] + (-depth*dt)*(1.0/dx*(u[1:,:]-u[:-1,:]) + 1.0/dy * (v[:,1:]-v[:,:-1]))

# Create an integrator function
f = SXFunction([p,uk,vk,hk],[u,v,h])
f.init()

# Mayer term
m = SXFunction([hk,hmeas],[sumAll((hk-hmeas)**2)])
m.init()

print "generated discrete dynamics"

# Allocate memory
u = copy.deepcopy(u0)
v = copy.deepcopy(v0)
h = copy.deepcopy(h0)

# Prepare plotting
if plot_progress:
  plt.figure(1)
  plt.clf()
  plt.grid(True)
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
    plt.ioff()
    plt.contour(X[::numboxes_per_plot],Y[::numboxes_per_plot],h[::numboxes_per_plot])
    plt.axis([0,poolwidth,0,poolwidth])
    plt.draw()
    time.sleep(0.1)
    
  # Integrate with Euler forward
  for j in range(num_eulersteps):
    u[1:-1,:] += dt*(-g/dx * (h[1:,:]-h[:-1,:]) - u[1:-1,:]*drag0)
    v[:,1:-1] += dt*(-g/dy * (h[:,1:]-h[:,:-1]) - v[:,1:-1]*drag0)
    h[:,:]    += (-depth0*dt)*(1.0/dx*(u[1:,:]-u[:-1,:]) + 1.0/dy * (v[:,1:]-v[:,:-1]))

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

# Solve with IPOPT
nlp_solver = IpoptSolver(ffcn)

# Set options
#nlp_solver.setOption("verbose",True)
nlp_solver.setOption("max_iter",10)
#nlp_solver.setOption("mu_init",1e-10)

# Initialize NLP solver
nlp_solver.init()

# Set initial condition and bounds
nlp_solver.setInput([2.0,0.01],NLP_X_INIT)
#nlp_solver.setInput([0.5,0.01],NLP_X_INIT)
nlp_solver.setInput([0,0],NLP_LBX)

# Solve single shooting problem
nlp_solver.evaluate()


