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
disp 'Testing sensitivity analysis in CasADi'

% All ODE and DAE integrators to be tested
DAE_integrators = {'idas','collocation'};
ODE_integrators = {'cvodes','rk', DAE_integrators{:}};

for ode=0:1
  if ode
    disp '******'
    disp 'Testing ODE example'
    Integrators = ODE_integrators;

    % Time
    t = SX.sym('t');

    % Parameter
    u = SX.sym('u');

    % Differential states
    s = SX.sym('s'); v = SX.sym('v'); m = SX.sym('m');
    x = [s;v;m];

    % Constants
    alpha = 0.05; % friction
    beta = 0.1;   % fuel consumption rate

    % Differential equation
    ode = [v;
           (u-alpha*v*v)/m;
           -beta*u*u];

    % Quadrature
    quad = v^3 + ((3-sin(t)) - u)^2;

    % DAE
    dae = struct('t', t, 'x', x, 'p', u, 'ode', ode, 'quad', quad);

    % Time length
    tf = 0.5;

    % Initial position
    x0 = [0;0;1];

    % Parameter
    u0 = 0.4;

  else
    disp '******'
    disp 'Testing DAE example'
    Integrators = DAE_integrators;

    % Differential state
    x = SX.sym('x');

    % Algebraic variable
    z = SX.sym('z');

    % Parameter
    u = SX.sym('u');

    % Differential equation
    ode = -x + 0.5*x*x + u + 0.5*z;

    % Algebraic constraint
    alg = z + exp(z) - 1.0 + x;

    % Quadrature
    quad = x*x + 3.0*u*u;

    % DAE
    dae = struct('x', x, 'z', z, 'p', u, 'ode', ode, 'alg', alg, 'quad', quad);

    % End time
    tf = 5;

    % Initial position
    x0 = 1;

    % Parameter
    u0 = 0.4;
  end

  % Integrator
  for integ=1:numel(Integrators)
    MyIntegrator = Integrators{integ};

    disp(sprintf('========'));
    disp(sprintf('Integrator: %s', MyIntegrator));
    disp(sprintf('========'));

    % Integrator options
    opts = struct('tf', tf);
    if strcmp(MyIntegrator,'collocation')
      opts.rootfinder = 'kinsol';
    end

    % Integrator
    I = casadi.integrator('I', MyIntegrator, dae, opts);

    % Integrate to get results
    res = I('x0', x0, 'p', u0);
    xf = full(res.xf);
    qf = full(res.qf);
    fprintf('%50s: xf=%s, qf=%s\n', 'Unperturbed solution', ...
            sprintf('%d ', xf), sprintf('%d ', qf));

    % Perturb solution to get a finite difference approximation
    h = 0.001;
    res = I('x0', x0, 'p', u0+h);
    fd_xf = (full(res.xf)-xf)/h;
    fd_qf = (full(res.qf)-qf)/h;
    fprintf('%50s: d(xf)/d(p)=%s, d(qf)/d(p)=%s\n', 'Finite differences', ...
            sprintf('%d ', fd_xf), sprintf('%d ', fd_qf));

    % Calculate one directional derivative, forward mode
    I_fwd = I.factory('I_fwd', {'x0', 'z0', 'p', 'fwd:p'}, {'fwd:xf', 'fwd:qf'});
    res = I_fwd('x0', x0, 'p', u0, 'fwd_p', 1);
    fwd_xf = full(res.fwd_xf);
    fwd_qf = full(res.fwd_qf);
    fprintf('%50s: d(xf)/d(p)=%s, d(qf)/d(p)=%s\n', 'Forward sensitivities', ...
            sprintf('%d ', fwd_xf), sprintf('%d ', fwd_qf));

    % Calculate one directional derivative, reverse mode
    I_adj = I.factory('I_adj', {'x0', 'z0', 'p', 'adj:qf'}, {'adj:x0', 'adj:p'});
    res = I_adj('x0', x0, 'p', u0, 'adj_qf', 1);
    adj_x0 = full(res.adj_x0);
    adj_p = full(res.adj_p);
    fprintf('%50s: d(qf)/d(x0)=%s, d(qf)/d(p)=%s\n', 'Adjoint sensitivities', ...
            sprintf('%d ', adj_x0), sprintf('%d ', adj_p));

    % Perturb adjoint solution to get a finite difference approximation of
    % the second order sensitivities
    res = I_adj('x0', x0, 'p', u0+h, 'adj_qf', 1);
    fd_adj_x0 = (full(res.adj_x0)-adj_x0)/h;
    fd_adj_p = (full(res.adj_p)-adj_p)/h;
    fprintf('%50s: d2(qf)/d(x0)d(p)=%s, d2(qf)/d(p)d(p)=%s\n', ...
            'FD of adjoint sensitivities', ...
            sprintf('%d ', fd_adj_x0), sprintf('%d ', fd_adj_p));

    % Forward over adjoint to get the second order sensitivities
    I_foa = I_adj.factory('I_foa', {'x0', 'z0', 'p', 'adj_qf', 'fwd:p'}, ...
                          {'fwd:adj_x0', 'fwd:adj_p'});
    res = I_foa('x0', x0, 'p', u0, 'adj_qf', 1, 'fwd_p', 1);
    fwd_adj_x0 = full(res.fwd_adj_x0);
    fwd_adj_p = full(res.fwd_adj_p);
    fprintf('%50s: d2(qf)/d(x0)d(p)=%s, d2(qf)/d(p)d(p)=%s\n', ...
            'Forward over adjoint sensitivities', ...
            sprintf('%d ', fd_adj_x0), sprintf('%d ', fd_adj_p));

    % Adjoint over adjoint to get the second order sensitivities
    I_aoa = I_adj.factory('I_aoa', {'x0', 'z0', 'p', 'adj_qf', 'adj:adj_p'}, ...
                         {'adj:x0', 'adj:p'});
    res = I_aoa('x0', x0, 'p', u0, 'adj_qf', 1, 'adj_adj_p', 1);
    adj_adj_x0 = full(res.adj_x0);
    adj_adj_p = full(res.adj_p);
    fprintf('%50s: d2(qf)/d(x0)d(p)=%s, d2(qf)/d(p)d(p)=%s\n', ...
            'Adjoint over adjoint sensitivities', ...
            sprintf('%d ', adj_adj_x0), sprintf('%d ', adj_adj_p));
  end
end
