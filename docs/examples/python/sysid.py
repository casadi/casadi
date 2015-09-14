from casadi import *

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
N = 10000  # Number of samples
fs = 610.1 # Sampling frequency [hz]

param_truth = DMatrix([5.625e-6,2.3e-4,1,4.69])
param_guess = DMatrix([5,2,1,5])
scale = vertcat([1e-6,1e-4,1,1])

############ MODELING #####################
y  = MX.sym('y')
dy = MX.sym('dy')
u  = MX.sym('u')

states = vertcat([y,dy]);
controls = u;

M = MX.sym("x")
c = MX.sym("c")
k = MX.sym("k")
k_NL = MX.sym("k_NL")

params = vertcat([M,c,k,k_NL])

rhs = vertcat([dy , (u-k_NL*y**3-k*y-c*dy)/M])

# Form an ode function
ode = MXFunction('ode',[states,controls,params],[rhs])

############ Creating a simulator ##########
N_steps_per_sample = 10
dt = 1/fs/N_steps_per_sample

# Build an integrator for this system: Runge Kutta 4 integrator
[k1] = ode([states,controls,params])
[k2] = ode([states+dt/2.0*k1,controls,params])
[k3] = ode([states+dt/2.0*k2,controls,params])
[k4] = ode([states+dt*k3,controls,params])

states_final = states+dt/6.0*(k1+2*k2+2*k3+k4)

# Create a function that simulates one step propagation in a sample
one_step = MXFunction('one_step',[states, controls, params],[states_final]);

X = states;

for i in range(N_steps_per_sample):
  [X] = one_step([X, controls, params])

# Create a function that simulates all step propagation on a sample
one_sample = MXFunction('one_sample',[states, controls, params], [X])

# speedup trick: expand into scalar operations
one_sample = one_sample.expand()

############ Simulating the system ##########

all_samples = one_sample.mapaccum("all_samples", N)

# Choose an excitation signal
numpy.random.seed(0)
u_data = DMatrix(0.1*numpy.random.random(N))

x0 = DMatrix([0,0])
[X_measured] = all_samples([x0, u_data, repmat(param_truth,1,N) ])

y_data = X_measured[0,:].T

# You may add some noise here
#y_data+= 0.001*numpy.random.random(N)
# When noise is absent, the fit will be perfect.

############ Create a Gauss-Newton solver ##########
def gauss_newton(e,nlp,V):
  gradF = nlp.gradient()
  jacG = nlp.jacobian("x","g")

  gradF.derivative(0, 1)

  J = jacobian(e,V)
  sigma = MX.sym("sigma")
  hessLag = MXFunction('H',hessLagIn(x=V,lam_f=sigma),hessLagOut(hess=sigma*mul(J.T,J)))

  # Code-generation for extra speedup
  nlp     = JitCompiler("clang",nlp)
  gradF   = JitCompiler("clang",gradF)
  jacG    = JitCompiler("clang",jacG)
  hessLag = JitCompiler("clang",hessLag)

  options = {}
  options["hess_lag"] = hessLag
  options["grad_f"] = gradF
  options["jac_g"] = jacG
  return NlpSolver("solver","ipopt", nlp, options)

############ Identifying the simulated system: single shooting strategy ##########

# Note, it is in general a good idea to scale your decision variables such
# that they are in the order of ~0.1..100
[X_symbolic] = all_samples([x0, u_data, repmat(params*scale,1,N) ])

e = y_data-X_symbolic[0,:].T;

nlp = MXFunction("nlp", nlpIn(x=params), nlpOut(f=0.5*inner_prod(e,e)))

solver = gauss_newton(e,nlp, params)

sol = solver(x0=param_guess)

print sol["x"]*scale

assert(max(fabs(sol["x"]*scale-param_truth))<1e-8)

############ Identifying the simulated system: multiple shooting strategy ##########

# All states become decision variables
X = MX.sym("X", 2, N)

[Xn] = one_sample.map([X, u_data.T, repmat(params*scale,1,N)])

gaps = Xn[:,:-1]-X[:,1:]

e = y_data-Xn[0,:].T;

V = veccat([params, X])

nlp = MXFunction("nlp", nlpIn(x=V), nlpOut(f=0.5*inner_prod(e,e),g=gaps))

# Multipleshooting allows for careful initialization
yd = np.diff(y_data,axis=0)*fs
X_guess = horzcat([ y_data , vertcat([yd,yd[-1]]) ]).T;

x0 = veccat([param_guess,X_guess])

solver = gauss_newton(e,nlp, V)

sol = solver(x0=x0,lbg=0,ubg=0)

print sol["x"][:4]*scale

assert(max(fabs(sol["x"][:4]*scale-param_truth))<1e-8)
