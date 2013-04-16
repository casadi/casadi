from casadi import *
from casadi.tools import *

from pylab import *

# System states
states = struct_ssym(["x","y","dx","dy"])
x,y,dx,dy = states[...]
    
# System controls
controls = struct_ssym(["u","v"])
u,v = controls[...]

# System parameters
parameters = struct_ssym(["k","c","beta"])
k,c,beta = parameters[...]

# Provide some numerical values
parameters_ = parameters()
parameters_["k"] = 10
parameters_["beta"] = 1
parameters_["c"] = 1

vel = vertcat([dx,dy])
p = vertcat([x,y])
q = vertcat([u,v])

# System dynamics
F = -k*(p-q) - beta*v*sqrt(sumAll(vel**2)+c**2)

# System right hand side
rhs = struct_SX(states)
rhs["x"]  = dx
rhs["y"]  = dy
rhs["dx"] = F[0]
rhs["dy"] = F[1]

f = SXFunction(controldaeIn(x=states,p=parameters,u=controls),daeOut(ode=rhs))
f.init()

N = 100

tgrid = linspace(0,10.0,N)
csim = ControlSimulator(f,tgrid)
csim.setOption("integrator",CVodesIntegrator)
csim.init()

x0 = states(0)

control_ = controls.repeated(csim.input("u"))
control_[0,"u"] = 1     # Kick the system with u=1 at the start
control_[N/2,"v"] = 2   # Kick the system with v=2 at half the simulation time

csim.setInput(x0,"x0")
csim.setInput(parameters_,"p")
csim.evaluate()

output = states.repeated(csim.output())

plot(tgrid,output[vertcat,:,"x"])
plot(tgrid,output[vertcat,:,"y"])

show()

