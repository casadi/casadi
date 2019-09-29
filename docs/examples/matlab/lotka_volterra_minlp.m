%
%     This file is part of CasADi.
%
%     CasADi -- A symbolic framework for dynamic optimization.
%     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
%                             K.U. Leuven. All rights reserved.
%     Copyright (C) 2011-2014 Greg Horn
%
%     CasADi is free software; you can redistribute it and/or
%     modify it under the terms of the GNU Lesser General Public
%     License as published by the Free Software Foundation; either
%     version 3 of the License, or (at your option) any later version.
%
%     CasADi is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     Lesser General Public License for more details.
%
%     You should have received a copy of the GNU Lesser General Public
%     License along with CasADi; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%

import casadi.*

% Lotka-Volerra problem from mintoc.de

T = 12; % Time horizon
N = 40; % number of control intervals

% Model parameters
c0 = 0.4;
c1 = 0.2;

% Initial condition for x
x0 = [0.5;0.7];

% Bounds on x
lbx = [0;0];
ubx = [2;2];

% Bounds on u
lbu = 0;
ubu = 1;

% Declare model variables
x = SX.sym('x', 2);
u = SX.sym('u');

% Model equations
xdot = [x(1) - x(1) * x(2) - c0 * x(1) * u;
              -x(2) + x(1) * x(2) - c1 * x(2) * u];

% Objective term
L = (x(1) - 1)^2 + (x(2) - 1)^2 + 1e-4*u^2;

% Formulate discrete time dynamics
if false
   % CVODES from the SUNDIALS suite
   dae = struct('x',x, 'p',u, 'ode',xdot, 'quad',L);
   opts = struct('tf',T/N);
   F = integrator('F', 'cvodes', dae, opts);
else
   % Fixed step Runge-Kutta 4 integrator
   M = 4; % RK4 steps per interval
   DT = T/N/M;
   f = Function('f', {x, u}, {xdot, L});
   X0 = MX.sym('X0', 2);
   U = MX.sym('U');
   X = X0;
   Q = 0;
   for j=1:M
       [k1, k1_q] = f(X, U);
       [k2, k2_q] = f(X + DT/2 * k1, U);
       [k3, k3_q] = f(X + DT/2 * k2, U);
       [k4, k4_q] = f(X + DT * k3, U);
       X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
       Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
   end
   F = Function('F', {X0, U}, {X, Q},char('x0','p'),char('xf','qf'));

end
% Initial guess for u
u_start = zeros(1,N);

% Get a feasible trajectory as an initial guess
xk = x0;
x_start = [xk];
for k=1:N
    ret = F('x0',xk, 'p',u_start(k));
    xk = ret.xf;
    x_start = [x_start xk];
end

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
discrete = [];
J = 0;
g={};
lbg = [];
ubg = [];

% "Lift" initial conditions
X0 = MX.sym('X0', 2);
w = {w{:} X0};
lbw = [lbw; x0];
ubw = [ubw; x0];
w0 = [w0; x_start(:,1)];
discrete = [discrete; 0; 0];

% Formulate the NLP
Xk = X0;
for k=1:N
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)]);
    w   = {w{:} Uk};
    lbw = [lbw;lbu];
    ubw = [ubw;ubu];
    w0  = [w0;u_start(:,k)];
    discrete = [discrete;1];

    % Integrate till the end of the interval
    Fk = F('x0',Xk,'p',Uk);
    Xk_end = Fk.xf;
    J=J+Fk.qf;

    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], 2);
    w   = {w{:} Xk};
    lbw = [lbw;lbx];
    ubw = [ubw;ubx];
    w0  = [w0;x_start(:,k+1)];
    discrete = [discrete;0;0];

    % Add equality constraint
    g   = {g{:} Xk_end-Xk};
    lbg = [lbg; 0; 0];
    ubg = [ubg; 0; 0];
end

% Concatenate decision variables and constraint terms
w = vertcat(w{:});
g = vertcat(g{:});

% Create an NLP solver
nlp_prob = struct('f', J, 'x', w, 'g', g);
nlp_solver = nlpsol('nlp_solver', 'bonmin', nlp_prob, struct('discrete', discrete));
%nlp_solver = nlpsol('nlp_solver', 'knitro', nlp_prob, {"discrete": discrete});
%nlp_solver = nlpsol('nlp_solver', 'ipopt', nlp_prob); % Solve relaxed problem

% Plot the solution
tgrid = 0:T/N:T;
figure(1)
clf()
hold on

% Solve the NLP
sol = nlp_solver('x0',w0, 'lbx',lbw, 'ubx',ubw, 'lbg',lbg, 'ubg',ubg);
w_opt = full(sol.x);
lam_w_opt = sol.lam_x;
lam_g_opt = sol.lam_g;
x0_opt = w_opt(1:3:end);
x1_opt = w_opt(2:3:end);
u_opt = w_opt(3:3:end);
plot(tgrid, x0_opt, '--');
plot(tgrid, x1_opt, '-');
stairs(tgrid, [nan; u_opt], '-.');
xlabel('t');
legend('x0','x1','u');
grid('on');

