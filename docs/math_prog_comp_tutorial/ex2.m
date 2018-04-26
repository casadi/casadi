% MATLAB/Octave
import casadi.*

% Formulate the NLP
x=SX.sym('x');
y=SX.sym('y');
z=SX.sym('z');
f=x^2+100*z^2;
g=z+(1-x)^2-y;
P=struct('x',[x;y;z],...
         'f',f,'g',g);

% Create solver instance
F=nlpsol('F','ipopt',P);

% Solve the problem
r=F('x0',[2.5 3.0 0.75],...
    'ubg',0,'lbg',0);
disp(r.x)
