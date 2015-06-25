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
    ode = [v,
           (u-alpha*v*v)/m,
           -beta*u*u];
      
    % Quadrature
    quad = v^3 + ((3-sin(t)) - u)^2;

    % DAE callback function
    ffcn = SXFunction('ffcn', daeIn('t',t,'x',x,'p',u),daeOut('ode',ode,'quad',quad));

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

    % DAE callback function
    ffcn = SXFunction('ffcn', daeIn('x',x,'z',z,'p',u),daeOut('ode',ode,'alg',alg,'quad',quad));
    
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
    disp(opts)

    % Integrator
    %I = Integrator('myintegrator', MyIntegrator, ffcn); % segfault
  end
end
