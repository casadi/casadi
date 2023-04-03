import pylab as plt
#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
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
# -*- coding: utf-8 -*-
from casadi import *

# Pendulum: point mass on massless rod 
# x'' = T*x, y'' = T*y-g, x**2+y**2-L**2
#
#
#     .
#  L /
#   /
#  /
# o (x,y) 
#
x = SX.sym("x") # horizontal position of point
y = SX.sym("y") # vertical position of point

dx = SX.sym("dx")
dy = SX.sym("dy")

u = SX.sym("u") # helper states to bring second order system to first order
v = SX.sym("v")

du = SX.sym("du")
dv = SX.sym("dv")

# Force
T = SX.sym("T")


L = 1
g = 9.81

alg = vertcat(dx-u,dy-v,du-T*x,dv-T*y+g,x**2+y**2-L**2)

# Implicit differential states
x_impl  = vertcat(x,y,u,v)
# Derivatives of implict differential states
dx_impl = vertcat(dx,dy,du,dv)

z = T

dae = {"x_impl": x_impl, "dx_impl": dx_impl, "z": z, "alg": alg}

# Perform index reduction
(dae_reduced,stats) = dae_reduce_index(dae)
print("We had an index 3 system (reduced to index 1 now)")
print(stats)
print("Here's the reduced DAE:")
print(dae_reduced)
print("Notice how the algebraic equation now contains a second derivative of x**2+y**2-L**2")
print("Also note how x**2+y**2-L**2 and its first derivate are absent from the algebraic equations.")
print("Both are kept in 'I' though: the resulting DAE invariants.")

# The DAE is not yet in a form that CasADi integrator can handle (semi-explicit).
# Let's convert it, ad obtain some adaptors
[dae_se, state_to_orig, phi] = dae_map_semi_expl(dae, dae_reduced)

grid = list(np.linspace(0,5,1000))
tf = DM(grid).T
intg = integrator("intg","idas",dae_se,0,grid)

# Before we can integrate, we need a consistant initial guess
# Let's say we start at y=-0.5, dx=-0.1;
# consistency equations should automatically identify everything else

# Encode the desire to only enforce y and dx
init_strength = {}
free  = 0 # default
force = -1
init_strength["x_impl"] = vertcat(free,force,free,free)
init_strength["dx_impl"] = vertcat(force,free,free,free)
init_strength["z"] = vertcat(free)

# Obtain a generating Function for initial guesses
init_gen = dae_init_gen(dae, dae_reduced, "ipopt", init_strength)

init = {}
blank = 0.0 # Will not be enforced but still used as initial guess
# suggest to initialize with the left-hanging solution by having x=-1 as initial guess 
init["x_impl"]  = vertcat(-1.0,-0.5,blank,blank)  
init["dx_impl"] = vertcat(-0.1,blank,blank,blank)
init["z"]  = blank

xz0 = init_gen(**init)
print("We have a consistent initial guess now to feed to the integrator:")
print(xz0)

print("Look, we found that force in the pendulum at t=0 equals:")
print(state_to_orig(xf=xz0['x0'],zf=xz0['z0'])["z"])

print("A consistent initial guess should make the invariants zero (up to solver precision):")
print(phi(x=xz0['x0'],z=xz0['z0']))

# Integrate and get resultant xf,zf on grid
sol = intg(**xz0)

print(sol['xf'].shape)

# Solution projected into original DAE space
sol_orig = state_to_orig(xf=sol["xf"],zf=sol["zf"])

# We can see the large-angle pendulum motion play out well in the u state
plt.plot(tf.T,sol_orig["x_impl"].T)
plt.grid(True)
plt.legend([e.name() for e in vertsplit(x_impl)])
plt.xlabel("Time [s]")
plt.title("Boundary value problem solution trajectory")

# A perfect integrator will perfectly preserve the values of the invariants over time
# Integrator errors make the invariants drift in practice
# This is not a detail; the pendulum length is in fact growing!
error = phi(x=sol["xf"],z=sol["zf"])["I"]

plt.figure()
plt.plot(tf.T,error.T)
plt.grid(True)
plt.xlabel("Time [s]")
plt.title("Evolution trajectory of invariants")

# There are some techniques to avoid the drifting of invariants associated with index reduction
#  - Method of Dummy Derviatives (not implemented)
#    hard implementation details: may need to switch between different choices online
#     -> integration in parts triggered by zero-crossing events
#  - Stop the integration every once in a while and project back (not implemented)
#  - Baumgarte stabilization (implemented): build into the equations a form of continuous feedback that brings back deviations in invariants back into the origin


# Demonstrate Baumgarte stabilization with a pole of -1.
# Drift is absent now.
(dae_reduced,stats) = dae_reduce_index(dae, {"baumgarte_pole": -1})
[dae_se, state_to_orig, phi] = dae_map_semi_expl(dae, dae_reduced)
intg = integrator("intg","idas",dae_se,0,grid)
init_gen = dae_init_gen(dae, dae_reduced, "ipopt", init_strength)
xz0 = init_gen(**init)
sol = intg(**xz0)
error = phi(x=sol["xf"],z=sol["zf"])["I"]
plt.figure()
plt.plot(tf.T,error.T)
plt.grid(True)
plt.xlabel("Time [s]")
plt.title("Boundary value problem solution trajectory with Baumgarte pole=-1")



# Sweep across different integrator precisions and Baumgarte poles.
# The choice of pole is not really straightforward
plt.figure()

poles = [-0.1,-1,-10]
precisions = [0.1,1,10]

abstol_default = 1e-8
reltol_default = 1e-6

plt.subplot(len(precisions),len(poles),1)

for i,pole in enumerate(poles):
  for j,precision in enumerate(precisions):
    (dae_reduced,stats) = dae_reduce_index(dae, {"baumgarte_pole": pole})
    [dae_se, state_to_orig, phi] = dae_map_semi_expl(dae, dae_reduced)
    intg = integrator("intg","idas",dae_se,0,grid,{"abstol":abstol_default*precision,"reltol":reltol_default*precision})
    init_gen = dae_init_gen(dae, dae_reduced, "ipopt", init_strength, {"ipopt.print_level":0,"print_time": False})
    xz0 = init_gen(**init)
    sol = intg(**xz0)
    nfevals = intg.stats()["nfevals"]
    error = phi(x=sol["xf"],z=sol["zf"])["I"]
    plt.subplot(len(precisions),len(poles),j*len(poles)+i+1)
    plt.plot(tf.T,error.T)
    plt.grid(True)
    plt.xlabel("Time [s]")
    plt.title("pole=%0.1f, abstol=%0.0e => nfevals=%d" % (pole,abstol_default*precision,nfevals))

plt.show()
