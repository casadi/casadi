from casadi import *
import matplotlib.pyplot as plt

# Declare variables (use simple, efficient DAG)
t = SX("t") # time
x=SX("x"); y=SX("y"); u=SX("u"); L=SX("cost")

# ODE right hand side function
f = [(1 - y*y)*x - y + u, x, x*x + y*y + u*u]
rhs = SXFunction([[t],[x,y,L],[u]],[f])

# Create an integrator (CVodes)
I = CVodesIntegrator(rhs)
I.setOption("abstol",1e-10) # abs. tolerance
I.setOption("reltol",1e-10) # rel. tolerance
I.setOption("steps_per_checkpoint",100)
I.setOption("stop_at_end",True)
I.init()

# All controls (use complex, general DAG)
NU = 20; U = MX("U",NU)

# The initial state (x=0, y=1, L=0)
X  = MX([0,1,0])

# Time horizon
T0 = MX(0);  TF = MX(20.0/NU)

# State derivative (not used)
XP = MX()

# Algebraic state (not used)
Z = MX()

# Build up a graph of integrator calls
for k in range(NU):
  [X,XP,Z] = I.call([T0,TF,X,U[k],XP,Z])

# Objective function: L(T)
obj = X[2]
F = MXFunction([U],[obj])

# Terminal constraints: 0<=[x(T);y(T)]<=0
eq = X[0:2]
G = MXFunction([U],[eq])

# Use the lifted newton method, transforming the single shooting problem to a multiple shooting problem 
lifted_newton = False

if lifted_newton:
  # Create a liftopt solver instance
  solver = LiftoptSolver(F,G) # TODO: debug! 
  
  # Set some options
  solver.setOption("optimizer","sqp")
  solver.setOption("lifted",True)
  
else:
  # Single shooting
  solver = IpoptSolver(F,G)
  solver.setOption("tol",1e-5)
  solver.setOption("hessian_approximation", \
                  "limited-memory")
  solver.setOption("max_iter",1000)
  solver.setOption("linear_solver","ma57")
  #solver.setOption("verbose",True)
  
# Initialize the NLP solver
solver.init()

# Set bounds and initial guess
solver.setInput(NU*[-0.75], NLP_LBX)
solver.setInput(NU*[1.0],NLP_UBX)
solver.setInput(NU*[0.0],NLP_X_INIT)
solver.setInput([0,0],NLP_LBG)
solver.setInput([0,0],NLP_UBG)

# Solve the problem
solver.solve()

# Plot the results
plt.plot(solver.output(NLP_X_OPT))
plt.show()
