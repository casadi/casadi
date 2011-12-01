from casadi import *
from casadi import tools
from matplotlib import pylab as plt

# Example 5.2 in Albersmeyer paper
# Solve F(u) = u**16 - 2 == 0

# Original variables
u = ssym("u")

# Lifted variables
x = ssym("x",4)

# Algorithm
xdef = SXMatrix.zeros(4)
xdef[0] = u**2
xdef[1] = x[0]**2
xdef[2] = x[1]**2
xdef[3] = x[2]**2
F = x[3] - 2

# Residual function G
G = SXFunction([u,x],[xdef-x,F])
G.init()

# Difference vector d
d = ssym("d",4)

# Substitute in the new definition of x
zdef = xdef-d
z = substituteInPlace(x,zdef)

# Eliminate x from F
F = substitute(F,x,z)

# Modified function Z
Z = SXFunction([u,d],[z,F])
Z.init()

# Matrix A and B in lifted Newton
A = Z.jac(0,0)
B = Z.jac(0,1) # calculate together!
AB  = SXFunction([u,d],[A,B])
AB.init()

# Variables
uk = 0.8*DMatrix.ones(1)
dk = DMatrix.zeros(4)
xk = DMatrix.nan(4)
Fk = DMatrix.nan()

# Stopping tolerance
TOL = 1e-10

# Initialize x0 by function evaluation
Z.setInput(uk,0)
Z.setInput(dk,1)
Z.evaluate()
Z.getOutput(xk,0)
Z.getOutput(Fk,1)

# Maximum number of iterations
max_iter = 100

# Print header
print " %4s" % "iter", " %20s" % "norm_Fk", " %20s" % "norm_dk"

# Iterate
k = 0
while True:
  
  # Get Ak and Bk
  AB.setInput(uk,0)
  AB.setInput(dk,1)
  AB.evaluate()
  Ak = AB.output(0)
  Bk = AB.output(1)
  
  # Get ak and bk
  Z.setInput(uk,0)
  Z.setInput(dk,1)
  Z.setFwdSeed(0,0)
  Z.setFwdSeed(dk,1)
  Z.evaluate(1,0)
  Z.getOutput(xk,0)
  Z.getOutput(Fk,1)
  ak = -Z.fwdSens(0)
  bk = Fk-Z.fwdSens(1)

  # Solve the condensed Newton system
  du = -solve(Bk,bk)
  
  # Perform the Newton step
  xk = xk + ak + mul(Ak,du)
  uk = uk + du
  
  # Call algorithm 2 to obtain new dk and Fk
  G.setInput(uk,0)
  G.setInput(xk,1)
  G.evaluate()
  G.getOutput(dk,0)
  G.getOutput(Fk,1)
  
  # Get error
  norm_Fk = float(norm_2(Fk))
  norm_dk = float(norm_2(dk))
  
  # Print
  print " %4d" % k, " %20e" % norm_Fk, " %20e" % norm_dk
  
  # Check if stopping criteria achieved
  if norm_Fk + norm_dk < TOL:
    print "Convergens achieved!"
    break
  
  # Increase iteration count
  k = k+1
  
  # Check if number of iterations have been reached
  if k >= max_iter:
    print "Maximum number of iterations (", max_iter, ") reached"
    break




#print A

#plt.spy(A.sparsity())
#plt.show()

#Z.init()
#print Z

