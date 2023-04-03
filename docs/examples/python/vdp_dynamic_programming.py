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

# End time
T = 10.

# Number of control intervals
N = 20

# Number of Runge-Kutta 4 steps per interval and step size
NK = 20
DT = T/(N*NK)

# Number of discrete control values
NU = 101

# Number of discrete state values
NX = 101

# System dynamics, can be called with matricex
def f(x1,x2,u):
  x1_dot = (1 - x2*x2)*x1 - x2 + u
  x2_dot = x1
  q_dot  = x1*x1 + x2*x2 + u*u
  return (x1_dot, x2_dot, q_dot)

# Control enumeration
U  = linspace(-1,1,NU)

# State space enumeration
x1 = linspace(-1,1,NX)
x2 = linspace(-1,1,NX)
X1,X2 = meshgrid(x1,x2)

# For each control action and state, precalculate next state and stage cost
stage_J = []
next_x1 = []
next_x2 = []
for u in U:
  # Take number of integration steps
  X1_k = copy(X1)
  X2_k = copy(X2)
  Q_k = zeros(X1.shape)
  for k in range(NK):
    # RK4 integration for x1, x2 and q
    k1_x1, k1_x2, k1_q = f(X1_k,                X2_k,                u)
    k2_x1, k2_x2, k2_q = f(X1_k + DT/2 * k1_x1, X2_k + DT/2 * k1_x2, u)
    k3_x1, k3_x2, k3_q = f(X1_k + DT/2 * k2_x1, X2_k + DT/2 * k2_x2, u)
    k4_x1, k4_x2, k4_q = f(X1_k + DT   * k3_x1, X2_k + DT   * k3_x2, u)
    X1_k += DT/6*(k1_x1 + 2*k2_x1 + 2*k3_x1 + k4_x1)
    X2_k += DT/6*(k1_x2 + 2*k2_x2 + 2*k3_x2 + k4_x2)
    Q_k  += DT/6*(k1_q  + 2*k2_q  + 2*k3_q  + k4_q )

  # Find out which state comes next (index)
  X1_k = matrix.round((X1_k+1)/2*(NX-1)).astype(int)
  X2_k = matrix.round((X2_k+1)/2*(NX-1)).astype(int)

  # Infinite cost if state gets out-of-bounds
  I = X1_k  <  0; Q_k[I]=inf; X1_k[I]=0
  I = X2_k  <  0; Q_k[I]=inf; X2_k[I]=0
  I = X1_k >= NX; Q_k[I]=inf; X1_k[I]=0
  I = X2_k >= NX; Q_k[I]=inf; X2_k[I]=0

  # Save the stage cost and next state
  next_x1.append(X1_k)
  next_x2.append(X2_k)
  stage_J.append(Q_k)

# Calculate cost-to-go (no end cost) and optimal control
J = zeros(X1.shape)
U_opt = []
for k in reversed(list(range(N))):
  # Cost to go for the previous step, optimal control action
  J_prev = inf*ones(X1.shape)
  u_prev = -ones(X1.shape,dtype=int)

  # Test all control actions
  for uind in range(NU):
    J_prev_test = J[next_x2[uind],next_x1[uind]]+stage_J[uind]
    better = J_prev_test<J_prev
    u_prev[better] = uind
    J_prev[better] = J_prev_test[better]

  # Update cost-to-go and save optimal control
  J = J_prev
  U_opt.append(u_prev)

# Reorder U_opt by stage
U_opt.reverse()

# Find optimal control starting at x1=0, x2=1
i1 = NX//2
i2 = NX-1
u_opt = []
x1_opt = [x1[i1]]
x2_opt = [x2[i2]]
cost = 0
for k in range(N):
  # Get the optimal control and go to next step
  u_ind = U_opt[k][i2,i1]
  cost += stage_J[u_ind][i2,i1]
  i1, i2 = next_x1[u_ind][i2,i1], next_x2[u_ind][i2,i1]

  # Save the trajectories
  u_opt.append(U[u_ind])
  x1_opt.append(x1[i1])
  x2_opt.append(x2[i2])

# Optimal cost
print("Minimal cost: ", cost)
assert abs(cost-J[NX-1,NX//2])<1e-8 # Consistency check

# Plot
figure(1)
clf()

# Plot optimal cost-to-go
subplot(121)
contourf(X1,X2,J)
colorbar()
xlabel('x1')
ylabel('x2')
title('Cost-to-go')

subplot(122)
plot(linspace(0,T,N+1),x1_opt,'--')
plot(linspace(0,T,N+1),x2_opt,'-.')
step(linspace(0,T,N),u_opt,'-')
plt.title("Dynamic programming solution")
plt.xlabel('time')
plt.legend(['x1 trajectory','x2 trajectory','u trajectory'])
grid(True)
show()
