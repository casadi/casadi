from casadi import *
from numpy import *
from matplotlib.pylab import plt
import time

# Physical parameters
g = 9.81
poolwidth = 0.2
drag = 4.0 # => u(0)
depth = 0.01 # => u(1)
sprad = 0.03
spheight = 0.01
endtime = 1.0

# Discretization
numboxes = 200
num_eulersteps = 100
num_measurements = 100

# Plotting
numboxes_per_plot = 1

# Discretization
ntimesteps = num_eulersteps*num_measurements
dt = endtime/ntimesteps
dx = poolwidth/numboxes
dy = poolwidth/numboxes
x = linspace(0,poolwidth,numboxes)
y = linspace(0,poolwidth,numboxes)
[X,Y] = meshgrid(x,y)

# Allocate memory
u = zeros((numboxes+1, numboxes))
v = zeros((numboxes  , numboxes+1))
h = zeros((numboxes  , numboxes))

# Initial conditions
spdist = sqrt( (X-0.04)**2 + (Y-0.04)**2)
I = spdist<sprad/3.0
h[I] = spheight * cos(3.0*pi*spdist[I]/(2.0*sprad))

# Clear figure
plt.figure(1)
plt.clf()
plt.grid(True)
plt.ion()
plt.hold(False)
plt.draw()
plt.show()

# Loop over all the measurements
for k in range(num_measurements):
  # Visualize the pool
  plt.ioff()
  plt.contour(X[::numboxes_per_plot],Y[::numboxes_per_plot],h[::numboxes_per_plot])
  plt.axis([0,poolwidth,0,poolwidth])
  plt.draw()
  time.sleep(0.1)
  
  # Integrate with Euler forward
  for j in range(num_eulersteps):
    u[1:-1,:] = u[1:-1,:] + dt*(-g/dx * (h[1:,:]-h[:-1,:]) - u[1:-1,:]*drag)
    v[:,1:-1] = v[:,1:-1] + dt*(-g/dy * (h[:,1:]-h[:,:-1]) - v[:,1:-1]*drag)
    h[:,:] = h[:,:] + (-depth*dt)*(1.0/dx*(u[1:,:]-u[:-1,:]) + 1.0/dy * (v[:,1:]-v[:,:-1]))

    