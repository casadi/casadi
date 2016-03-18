import casadi.*

% In this example, we fit a nonlinear model to measurements
%
% This example uses more advanced constructs than the vdp* examples:
% Since the number of control intervals is potentially very large here,
% we use memory-efficient Map and MapAccum, in combination with
% codegeneration.
%
% We will be working with a 2-norm objective:
% || y_measured - y_simulated ||_2^2
%
% This form is well-suited for the Gauss-Newton Hessian approximation.

%%%%%%%%%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%
N = 10000;  % Number of samples
fs = 610.1; % Sampling frequency [hz]

param_truth = [5.625e-6;2.3e-4;1;4.69];
param_guess = [5;2;1;5];
scale = [1e-6;1e-4;1;1];

%%%%%%%%%%%% MODELING %%%%%%%%%%%%%%%%%%%%%
y  = MX.sym('y');
dy = MX.sym('dy');
u  = MX.sym('u');

states = [y;dy];
controls = u;

M = MX.sym('x');
c = MX.sym('c');
k = MX.sym('k');
k_NL = MX.sym('k_NL');

params = [M;c;k;k_NL];

rhs = [dy; (u-k_NL*y.^3-k*y-c*dy)/M];

% Form an ode function
ode = Function('ode',{states,controls,params},{rhs});

%%%%%%%%%%%% Creating a simulator %%%%%%%%%%
N_steps_per_sample = 10;
dt = 1/fs/N_steps_per_sample;

% Build an integrator for this system: Runge Kutta 4 integrator
k1 = ode(states,controls,params);
k2 = ode(states+dt/2.0*k1,controls,params);
k3 = ode(states+dt/2.0*k2,controls,params);
k4 = ode(states+dt*k3,controls,params);

states_final = states+dt/6.0*(k1+2*k2+2*k3+k4);

% Create a function that simulates one step propagation in a sample
one_step = Function('one_step',{states, controls, params},{states_final});

X = states;
for i=1:N_steps_per_sample
    X = one_step(X, controls, params);
end

% Create a function that simulates all step propagation on a sample
one_sample = Function('one_sample',{states, controls, params}, {X});

% speedup trick: expand into scalar operations
one_sample = one_sample.expand();

%%%%%%%%%%%% Simulating the system %%%%%%%%%%

all_samples = one_sample.mapaccum('all_samples', N);

% Choose an excitation signal
u_data = 0.1*rand(N,1);

x0 = DM([0,0]);
X_measured = all_samples(x0, u_data, repmat(param_truth,1,N));

y_data = X_measured(1,:)';

%%%%%%%%%%%% Identifying the simulated system: single shooting strategy %%%%%%%%%%

% Note, it is in general a good idea to scale your decision variables such
% that they are in the order of ~0.1..100
X_symbolic = all_samples(x0, u_data, repmat(params.*scale,1,N));

e = y_data-X_symbolic(01,:)';

nlp = struct('x', params, 'f', 0.5*dot(e,e));
solver = nlpsol('solver','ipopt', nlp);

sol = solver('x0',param_guess);

sol.x.*scale

assert(max(full(abs(sol.x.*scale-param_truth)))<1e-8)
