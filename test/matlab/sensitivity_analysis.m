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
DAE_integrators = {'idas','collocation','oldcollocation'};
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
    opts = struct;
    opts.tf = tf;
    if strcmp(MyIntegrator,'collocation') | strcmp(MyIntegrator,'oldcollocation')
      opts.implicit_solver = 'kinsol';
      opts.implicit_solver_options = struct('linear_solver', 'csparse');
      if strcmp(MyIntegrator,'oldcollocation')
        opts.expand_f=true;
      end
    end

    % Integrator
    I = Function.ivpsol('myintegrator', MyIntegrator, dae, opts);

    % Integrate to get results
    arg = struct;
    arg.x0 = x0;
    arg.p = u0;
    res = I(arg);
    xf = full(res.xf);
    qf = full(res.qf);
    fprintf('Unperturbed solution: xf=%d, qf=%d\n', xf, qf);

    % Perturb solution to get a finite difference approximation
    h = 0.001;
    arg.p = u0+h;
    res = I(arg);
    xf_pert = full(res.xf);
    qf_pert = full(res.qf);
    fprintf('Finite differences: d(xf)/d(p)=%d, d(qf)/d(p)=%d\n',...
        (xf_pert-xf)/h, (qf_pert-qf)/h);

    % Calculate once directional derivative, forward mode
    I_fwd = I.forward(1);
    arg = struct('der_x0',x0,'der_p',u0,'der_xf',xf,'der_qf',qf);
    arg.fwd0_x0 = 0;
    arg.fwd0_p = 1;
    res = I_fwd(arg);
    fwd_xf = full(res.fwd0_xf);
    fwd_qf = full(res.fwd0_qf);
    fprintf('Forward sensitivities: d(xf)/d(p)=%d, d(qf)/d(p)=%d\n', ...
        fwd_xf, fwd_qf);

    % Calculate one directional derivative, reverse mode
    I_adj = I.reverse(1);
    arg = struct('der_x0',x0,'der_p',u0,'der_xf',xf,'der_qf',qf);
    arg.adj0_xf = 0;
    arg.adj0_qf = 1;
    res = I_adj(arg);
    adj_x0 = full(res.adj0_x0);
    adj_p = full(res.adj0_p);
    fprintf('Adjoint sensitivities: d(qf)/d(x0)=%d, d(qf)/d(p)=%d\n', ...
        adj_x0, adj_p);

    % Perturb adjoint solution to get a finite difference approximation of
    % the second order sensitivities
    arg.der_p = u0+h;
    res = I_adj(arg);
    adj_x0_pert = full(res.adj0_x0);
    adj_p_pert = full(res.adj0_p);
    fprintf('FD of adjoint sensitivities: d2(qf)/d(x0)d(p)=%d, d2(qf)/d(p)d(p)=%d\n', ...
        (adj_x0_pert-adj_x0)/h, (adj_p_pert-adj_p)/h);

    % Forward over adjoint to get the second order sensitivities
    I_foa = I_adj.forward(1);
    arg = struct('der_der_x0',x0,'der_der_p',u0,'der_der_xf',xf,'der_der_qf',qf);
    arg.der_adj0_x0 = adj_x0;
    arg.der_adj0_p = adj_p;
    arg.der_adj0_xf = 0;
    arg.der_adj0_qf = 1;
    arg.fwd0_der_p = 1.0;
    res = I_foa(arg);
    fwd_adj_x0 = full(res.fwd0_adj0_x0);
    fwd_adj_p = full(res.fwd0_adj0_p);
    fprintf('Forward over adjoint sensitivities: d2(qf)/d(x0)d(p)=%d, d2(qf)/d(p)d(p)=%d\n', ...
        fwd_adj_x0, fwd_adj_p);

    % Adjoint over adjoint to get the second order sensitivities
    I_aoa = I_adj.reverse(1);
    arg = struct('der_der_x0',x0,'der_der_p',u0,'der_der_xf',xf,'der_der_qf',qf);
    arg.der_adj0_x0 = adj_x0;
    arg.der_adj0_p = adj_p;
    arg.der_adj0_xf = 0;
    arg.der_adj0_qf = 1;
    arg.adj0_adj0_p = 1;
    res = I_aoa(arg);
    adj_adj_x0 = full(res.adj0_der_x0);
    adj_adj_p = full(res.adj0_der_p);
    fprintf('Adjoint over adjoint sensitivities: d2(qf)/d(x0)d(p)=%d, d2(qf)/d(p)d(p)=%d\n', ...
        adj_adj_x0, adj_adj_p);
  end
end

