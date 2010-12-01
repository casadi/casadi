function biegler_10_2

% optimal control problem
ocp = struct;

% Time 
ocp.t = sx_matrix('t');

% Differential states
Ls = sx_matrix('Ls');    % mean crystal size
Nc = sx_matrix('Nc');    % number of nuclei per liter of solvent
L  = sx_matrix('L');     % total length of crystals per liter of solvent
Ac = sx_matrix('Ac');    % total surface area of the crystals per liter of solvent
Vc = sx_matrix('Vc');    % total volume of the crysals per liter of solvent
Mc = sx_matrix('Mc');    % total mass of the crystals
Cc = sx_matrix('Cc');    % solute concentration
Tc = sx_matrix('Tc');    % cystillizer temperature

% State vector
ocp.x = [Ls;Nc;L;Ac;Vc;Mc;Cc;Tc];
ocp.nx = length(ocp.x);

% Bounds on the states
ocp.x_lb = inf(ocp.nx,1);
ocp.x_ub = inf(ocp.nx,1);

% Initial values
Ls_init = 0.0005;
Nc_init = 0;
L_init = 0;
Ac_init = 0;
Vc_init = 0;
Mc_init = 2.0;
Cc_init = 5.4;
Tc_init = 75;
ocp.x_init = [Ls_init;Nc_init;L_init;Ac_init;Vc_init;Mc_init;Cc_init;Tc_init];

% Control
Tj = sx_matrix('Tj'); % jacket temperature
Tj_lb = 10; Tj_ub = 60; Tj_init = 30;

% Control vector
ocp.u = Tj;
ocp.u_lb = Tj_lb;
ocp.u_ub = Tj_ub;
ocp.u_init = Tj_init;
ocp.nu = 1;

% Constants
Vs = 300; % volume of the solvent
W = 2025; % the total mass in the crysallizer
a = [-66.4309, 2.8604, -0.022579, 6.7117e-5];
b = [16.08852, -2.708263, 0.0670694, -3.5685e-4];
Kg = 0.00418;
Bn = 385;
Cp = 0.4;
Kc = 35;
Ke = 377;
eta1 = 1.1;
eta2 = 5.72;
Ls0 = 5e-4; % initial crystal size
L0 = 5e-5; % nucleate crystal size
Ws0 = 2; % weight of the seed crystals
rho = 1.58; % specific gravity of crystals
alpha = 0.2; % shape factor for area of crystals
beta = 1.2; % shape factor for volume of crystals

% Time horizon
ocp.tf = 16;

% Dependent variables
C_bar = 100*Cc/(1.35+Cc);
Tequ = a(1) + a(2)*C_bar + a(3)*C_bar*C_bar + a(4)*C_bar*C_bar*C_bar; % equilibrium temperature
Ta = b(1) + b(2)*C_bar + b(3)*C_bar*C_bar + b(4)*C_bar*C_bar*C_bar; % lower bound of Tj
% degree of supercooling
% DeltaT = max(0,Tequ-Tc);   % Original formulation
DeltaT = max(1e-8,Tequ-Tc); % 'epsilon' to avoid divide by zero
%  DeltaT = log(exp(0)+exp(Tequ-Tc));   % Log sum exp

% Differential equations
Ls_dot = Kg*sqrt(Ls)*DeltaT^eta1;
Nc_dot = Bn*DeltaT^eta2;
L_dot = Nc*Ls_dot + L0*Nc_dot;
Ac_dot = 2*alpha*Nc*Ls_dot + L0*L0*Nc_dot;
Vc_dot = 3*beta*Ac*Ls_dot + L0*L0*L0*Nc_dot;
Mc_dot = 3*(Ws0/(Ls0*Ls0*Ls0))*Ls*Ls*Ls_dot + rho*Vs*Vc_dot;
Cc_dot = -1/Vs*Mc_dot;
Tc_dot = (Kc*Mc_dot - Ke*(Tc-Tj))/(W*Cp);
xdot = [Ls_dot; Nc_dot; L_dot; Ac_dot; Vc_dot; Mc_dot; Cc_dot; Tc_dot];

% Right hand side of the ODE
y = {ocp.t,ocp.x,ocp.u};
ocp.ffcn = sx_function(y,xdot);

% Objective function (meyer term)
ocp.mfcn = sx_function(y,-Ls);

% Nonlinear constraint function
ocp.cfcn = sx_function(y,Tj-Ta);

% Degree of interpolating polynomial
K = 3;
  
% Number of finite elements
N = 15;

% Radau collocation points
cp = 'radau';

% Size of the finite elements
h = ocp.tf/N;

% Roots of the lagrange polynomials
tau_roots = get_collocation_points(K,cp);

% Coefficients of the collocation and continuity equations
[C,D] = get_collocation_coeff(K, cp);

% Collocated times
T = cell(N,K+1);
for i=1:N
    for j=1:K+1
        T{i,j} = h*(i-1 + tau_roots(j));
    end
end

% Collocated states
X = collocate('X',ocp.nx,N,K);

% Collocated control (piecewice constant)
U = collocate('U',ocp.nu,N);

% State at end time
XF = sx_matrix('XF',ocp.nx);

% All variables with bounds and initial guess
vars = sx_matrix();
vars_lb = [];
vars_ub = [];
vars_sol = [];

% Loop over the finite elements
for i=1:N
    % collocated controls 
    vars = [vars;U{i}];
    vars_lb = [vars_lb; ocp.u_lb];
    vars_ub = [vars_ub; ocp.u_ub];
    vars_sol = [vars_sol; ocp.u_init];

    % Collocated states
    for j = 1:K+1
        % Add to list of NLP variables
        vars = [vars;X{i,j}];
        vars_sol = [vars_sol; ocp.x_init];
        if (i==1 && j==1)
            % Initial constraints
            vars_lb = [vars_lb; ocp.x_init];
            vars_ub = [vars_ub; ocp.x_init];
        else
            % Variable bounds
            vars_lb = [vars_lb; ocp.x_lb];
            vars_ub = [vars_ub; ocp.x_ub];
        end
    end
end

% Add states at end time
vars = [vars;XF];
vars_sol = [vars_sol; ocp.x_init];
vars_lb = [vars_lb; ocp.x_lb];
vars_ub = [vars_ub; ocp.x_ub];
vars_sol = [vars_sol; ocp.x_init];

% Nonlinear constraint function for the NLP
g = sx_matrix();
lbg = [];
ubg = [];

for i=1:N
    for k=1:K
        % augmented state vector
        y_ik = {T{i,k}, X{i,k}, U{i}};

        % Add collocation equations to NLP
        rhs = sx_matrix(ocp.nx,1);
        for j=1:K
            rhs = rhs +  X{i,j}*C(j,k);
        end
        
        d = ocp.ffcn.evaluate_symbolic(y_ik);
        g = [g; h*d-rhs]; 
        lbg = [lbg; zeros(ocp.nx,1)];
        ubg = [ubg; zeros(ocp.nx,1)];

        % Add nonlinear constraints
        g = [g; ocp.cfcn.evaluate_symbolic(y_ik)];
        lbg = [lbg; 0];
        ubg = [ubg; inf];
    end
        
    % Add continuity equation to NLP
    rhs = sx_matrix(ocp.nx,1);
    for j=1:K+1
        rhs = rhs + D(j)*X{i,j};
    end
    if(i<N)
        g = [g; X{i+1,1} - rhs];
    else
        g = [g; XF - rhs];
    end
    lbg = [lbg; zeros(ocp.nx,1)];
    ubg = [ubg; zeros(ocp.nx,1)];
end

gfcn = sx_function(vars,g);

% Objective function of the NLP
y_f = {T{end,end},XF,U{end}};
f = ocp.mfcn.evaluate_symbolic(y_f);
ffcn = sx_function(vars, f);

% Jacobian of the constraints
jfcn = gfcn.jacobian();

% Hessian of the Lagrangian:
% Lagrange multipliers
lambda = sx_matrix('lambda',length(g));

% Objective function scaling
sigma = sx_matrix('sigma');

% Lagrangian function
lfcn = sx_function({vars,lambda,sigma}, sigma*f + lambda'*g);

% Hessian of the Lagrangian
hfcn = lfcn.hessian();

% ----
% SOLVE THE NLP
% ----

% Create a solver
solver = ipopt_solver(ffcn,gfcn,hfcn,jfcn);

% Set options
solver.set_option('abstol',1e-6);

% initialize the solver
solver.init();

% Initial condition
NLP_X_INIT = 0;
NLP_LBX = 1;
NLP_UBX = 2;
NLP_LBG = 3;
NLP_UBG = 4;
NLP_LAMBDA_INIT = 5;

NLP_X_OPT = 0;
NLP_COST = 1;
NLP_LAMBDA_OPT = 2;
NLP_LAMBDA_LBX = 3;
NLP_LAMBDA_UBX = 4;

solver.set_input(vars_sol,NLP_X_INIT);

% Bounds on x
solver.set_input(vars_lb,NLP_LBX);
 solver.set_input(vars_ub,NLP_UBX);

% Bounds on g
solver.set_input(lbg,NLP_LBG);
solver.set_input(ubg,NLP_UBG);

% Solve the problem
solver.evaluate();

% Print the optimal cost
cost = solver.get_output(NLP_COST)

% Get the solution
vars_sol = solver.get_output(NLP_X_OPT);


%   % Get the optimal solution
%   vector<double> Tj_opt(N);
%   vector<double> Ls_opt(N*(K+1));
%   vector<double> Nc_opt(N*(K+1));
%   vector<double> L_opt(N*(K+1));
%   vector<double> Ac_opt(N*(K+1));
%   vector<double> Vc_opt(N*(K+1));
%   vector<double> Mc_opt(N*(K+1));
%   vector<double> Cc_opt(N*(K+1));
%   vector<double> Tc_opt(N*(K+1));
%   vector<double> t_opt(N*(K+1));
%   int ind = 0; % index of nlp->x
%   for(int i=0; i<N; ++i){
%     Tj_opt[i] = vars_sol[ind++];
%     for(int j=0; j<=K; ++j){
%       int ij = (K+1)*i+j;
%       Ls_opt[ij] = vars_sol[ind++];
%       Nc_opt[ij] = vars_sol[ind++];
%       L_opt[ij] = vars_sol[ind++];
%       Ac_opt[ij] = vars_sol[ind++];
%       Vc_opt[ij] = vars_sol[ind++];
%       Mc_opt[ij] = vars_sol[ind++];
%       Cc_opt[ij] = vars_sol[ind++];
%       Tc_opt[ij] = vars_sol[ind++];
%       t_opt[ij] = T[i][j];
%     }
%   }
% 

  
end

function tau = get_collocation_points(K,cp)
switch(cp)
    case 'legendre'
        switch(K)
            case 1, tau = [0,0.500000];
            case 2, tau = [0,0.211325,0.788675];
            case 3, tau = [0,0.112702,0.500000,0.887298];
            case 4, tau = [0,0.069432,0.330009,0.669991,0.930568];
            case 5, tau = [0,0.046910,0.230765,0.500000,0.769235,0.953090];
            otherwise, error('wrong number of collocation points');
        end
    case 'radau'
        switch(K)
            case 1, tau = [0,1.000000];
            case 2, tau = [0,0.333333,1.000000];
            case 3, tau = [0,0.155051,0.644949,1.000000];
            case 4, tau = [0,0.088588,0.409467,0.787659,1.000000];
            case 5, tau = [0,0.057104,0.276843,0.583590,0.860240,1.000000];
            otherwise, error('wrong number of collocation points');
        end
    otherwise
        error('cp must be radau or legendre');    
end
end

function VAR = collocate(name,nvar,varargin)
if nargin==3
    N = varargin{1};
    VAR = cell(N,1);
    
  	% Loop over the finite elements
    for i=1:N
        VAR{i} = sx_matrix([name,'_',num2str(i)],nvar);
    end
else
    assert(nargin==4)
    N = varargin{1};
    K = varargin{2};
    VAR = cell(N,K+1);
    
    % Loop over the finite elements
    for i=1:N
        for j=1:K+1
            VAR{i,j} = sx_matrix([name,'_',num2str(i),'_',num2str(j)],nvar);
        end
    end
    
end
end

% Get the coefficeints for the collocation and continuity equations
function [C,D] = get_collocation_coeff(K, cp)
    D = nan(K+1,1);
    C = nan(K+1,K+1);
  
    % Collocation point
    tau = sx_matrix('tau');
  
    % Roots
    tau_root = get_collocation_points(K,cp);

    % Lagrange polynomials
    l = cell(K+1,1);
    for j=1:K+1
        L = sx_matrix(1);
        for k=1:K+1
            if(k ~= j)
                L = L * (tau-tau_root(k))/(tau_root(j)-tau_root(k));
                l{j} = sx_function(tau,L);
                l{j}.set_option('ad_order',1);
                l{j}.init();
            end
        end
    end
  
    % Get the coefficients of the continuity equation
    for j = 1:K+1
        l{j}.set_input(1);
        l{j}.evaluate();
        D(j) = l{j}.get_output();
    end

    % Get the coefficients of the collocation equation
    for j=1:K+1
        for k=1:K+1
            l{j}.set_input(tau_root(k));
            l{j}.evaluate();
      
            l{j}.set_input(1,0,1);
            l{j}.evaluate_fwd();
            C(j,k) = l{j}.get_output(0,1);
        end
    end
end