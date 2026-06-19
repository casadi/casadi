#
#     MIT No Attribution
#
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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
import casadi as ca
import numpy

# In this example, we fit a nonlinear model to measurements
#
# This example uses more advanced constructs than the vdp* examples:
# Since the number of control intervals is potentially very large here,
# we use memory-efficient Map and MapAccum, in combination with
# codegeneration.
#
# We will be working with a 2-norm objective:
# || y_measured - y_simulated ||_2^2
#
# This form is well-suited for the Gauss-Newton Hessian approximation.

############ SETTINGS #####################
N = 2000  # Number of samples
fs = 610.1 # Sampling frequency [hz]

param_truth = ca.DM([5.625e-6,2.3e-4,1,4.69])
param_guess = ca.DM([5,2,1,5])
scale = ca.vertcat(1e-6,1e-4,1,1)

############ MODELING #####################
y  = ca.MX.sym('y')
dy = ca.MX.sym('dy')
u  = ca.MX.sym('u')

states = ca.vertcat(y,dy)
controls = u

M = ca.MX.sym("M")
c = ca.MX.sym("c")
k = ca.MX.sym("k")
k_NL = ca.MX.sym("k_NL")

params = ca.vertcat(M,c,k,k_NL)

rhs = ca.vertcat(dy , (u-k_NL*y**3-k*y-c*dy)/M)

# Form an ode function
ode = ca.Function('ode',[states,controls,params],[rhs])

############ Creating a simulator ##########
N_steps_per_sample = 10
dt = 1/fs/N_steps_per_sample

# Build an integrator for this system: Runge Kutta 4 integrator
k1 = ode(states,controls,params)
k2 = ode(states+dt/2.0*k1,controls,params)
k3 = ode(states+dt/2.0*k2,controls,params)
k4 = ode(states+dt*k3,controls,params)

states_final = states+dt/6.0*(k1+2*k2+2*k3+k4)

# Create a function that simulates one step propagation in a sample
one_step = ca.Function('one_step',[states, controls, params],[states_final])

X = states

for i in range(N_steps_per_sample):
  X = one_step(X, controls, params)

# Create a function that simulates all step propagation on a sample
one_sample = ca.Function('one_sample',[states, controls, params], [X])

############ Simulating the system ##########
all_samples = one_sample.mapaccum("all_samples", N)

# Choose an excitation signal
numpy.random.seed(0)
u_data = ca.DM(0.1*numpy.random.random(N))

x0 = ca.DM([0,0])
X_measured = all_samples(x0, u_data, ca.repmat(param_truth,1,N))

y_data = X_measured[0,:].T

# You may add some noise here
#y_data+= 0.001*numpy.random.random(N)
# When noise is absent, the fit will be perfect.

# Use just-in-time compilation to speed up the evaluation
if ca.Importer.has_plugin('clang'):
  with_jit = True
  compiler = 'clang'
elif ca.Importer.has_plugin('shell'):
  with_jit = True
  print("WARNING: on Windows, JIT may require a special environment cfr https://github.com/casadi/casadi/wiki/FAQ:-how-to-perform-jit-for-function-evaluations-of-my-optimization-problem%3F")
  compiler = 'shell'
else:
  print("WARNING; running without jit. This may result in very slow evaluation times")
  with_jit = False
  compiler = ''

############ Create a Gauss-Newton solver ##########
def gauss_newton(e,nlp,V):
  J = ca.jacobian(e,V)
  H = ca.triu(J.T @  J)
  sigma = ca.MX.sym("sigma")
  hessLag = ca.Function('nlp_hess_l',{'x':V,'lam_f':sigma, 'hess_gamma_x_x':sigma*H},
                     ['x','p','lam_f','lam_g'], ['hess_gamma_x_x'],
                     dict(jit=with_jit, compiler=compiler))
  return ca.nlpsol("solver","ipopt", nlp, dict(hess_lag=hessLag, jit=with_jit, compiler=compiler))

############ Identifying the simulated system: single shooting strategy ##########

# Note, it is in general a good idea to scale your decision variables such
# that they are in the order of ~0.1..100
X_symbolic = all_samples(x0, u_data, ca.repmat(params*scale,1,N))

e = y_data-X_symbolic[0,:].T
nlp = {'x':params, 'f':0.5*ca.dot(e,e)}

solver = gauss_newton(e,nlp, params)

sol = solver(x0=param_guess)

print(sol["x"]*scale)

assert ca.norm_inf(sol["x"]*scale-param_truth)<1e-8

############ Identifying the simulated system: multiple shooting strategy ##########

# All states become decision variables
X = ca.MX.sym("X", 2, N)

Xn = one_sample.map(N, 'openmp')(X, u_data.T, ca.repmat(params*scale,1,N))

gaps = Xn[:,:-1]-X[:,1:]

e = y_data-Xn[0,:].T

V = ca.veccat(params, X)

nlp = {'x':V, 'f':0.5*ca.dot(e,e), 'g':ca.vec(gaps)}

# Multipleshooting allows for careful initialization
yd = numpy.diff(y_data,axis=0)*fs
X_guess = ca.horzcat(y_data , ca.vertcat(yd,yd[-1])).T

x0 = ca.veccat(param_guess,X_guess)

solver = gauss_newton(e,nlp, V)

sol = solver(x0=x0,lbg=0,ubg=0)

print(sol["x"][:4]*scale)

assert ca.norm_inf(sol["x"][:4]*scale-param_truth)<1e-8
