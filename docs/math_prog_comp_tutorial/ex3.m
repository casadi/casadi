% MATLAB/Octave
import casadi.*

% Formulate the ODE
x=SX.sym('x',2);
p=SX.sym('p');
z=1-x(2)^2;
f=[z*x(1)-x(2)+p;x(1)];
dae=struct('x',x,'p',p,...
           'ode',f);

% Create solver instance
op=struct('t0',0,'tf',1);
F=integrator('F',...
      'cvodes',dae,op);

% Solve the problem
r=F('x0',[0,1],'p',0.1);
disp(r.xf)

% Create Jacobian function
D=F.factory('D',...
 {'x0','p'},{'jac:xf:x0'});

% Solve the problem
r=D('x0',[0,1],'p',0.1);
disp(r.jac_xf_x0)
