from casadi import *
from casadi import tools
from matplotlib import pylab as plt

# Example 5.2 in Albersmeyer paper
# Solve F(u) = u**16 - 2 == 0

# Original variables
u = ssym("u")

# Algorithm
x1 = u**2
x2 = x1**2
x3 = x2**2
x4 = x3**2
F = x4 - 2

# Residual function
ffcn = SXFunction([u],[F])

# Lifting function
ifcn = SXFunction([u],[vertcat([x1,x2,x3,x4])])

# Initial guess for u
u_guess = 0.8


# Problem formulation ends
# Everything below should go into a lifted newton solver class

# Options
TOL = 1e-10     # Stopping tolerance
max_iter = 100  # Maximum number of iterations

# Extract the free variable and expressions for F and xdef
u = ffcn.inputSX()
f = ffcn.outputSX()
xdef = ifcn.outputSX()

# Lifted variables
x = ssym("x",xdef.size())

# Substitute in the lifted variables x into the expressions for xdef and F
ex = SXMatrixVector(1,f)
substituteInPlace(x, xdef, ex, True,False)
f = ex[0]

# Residual function G
G = SXFunction([u,x],[xdef-x,f])
G.init()

# Difference vector d
d = ssym("d",xdef.size())

# Substitute out the x from the zdef
zdef = xdef-d
z = substituteOut(x,zdef)

# Eliminate x from F
f = substitute(f,x,z)

# Modified function Z
Z = SXFunction([u,d],[z,f])
Z.init()

# Matrix A and B in lifted Newton
A = Z.jac(0,0)
B = Z.jac(0,1) # calculate together!
AB  = SXFunction([u,d],[A,B])
AB.init()

# Variables
uk = u_guess*DMatrix.ones(1)
dk = DMatrix.zeros(xdef.size())
xk = DMatrix.nan(xdef.size())
fk = DMatrix.nan()

# Initialize x0 by function evaluation
Z.setInput(uk,0)
Z.setInput(dk,1)
Z.evaluate()
Z.getOutput(xk,0)
Z.getOutput(fk,1)


# Print header
print " %4s" % "iter", " %20s" % "norm_fk", " %20s" % "norm_dk"

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
  Z.getOutput(fk,1)
  ak = -Z.fwdSens(0)
  bk = fk-Z.fwdSens(1)

  # Solve the condensed Newton system
  du = -solve(Bk,bk)
  
  # Perform the Newton step
  xk = xk + ak + mul(Ak,du)
  uk = uk + du
  
  # Call algorithm 2 to obtain new dk and fk
  G.setInput(uk,0)
  G.setInput(xk,1)
  G.evaluate()
  G.getOutput(dk,0)
  G.getOutput(fk,1)
  
  # Get error
  norm_fk = float(norm_2(fk))
  norm_dk = float(norm_2(dk))
  
  # Print
  print " %4d" % k, " %20e" % norm_fk, " %20e" % norm_dk
  
  # Check if stopping criteria achieved
  if norm_fk + norm_dk < TOL:
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

