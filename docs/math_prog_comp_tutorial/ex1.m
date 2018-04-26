% MATLAB/Octave
import casadi.*

% Symbolic representation
x=SX.sym('x');
y=SX.sym('y');
z=y-(1-x)^2;
f=x^2+100*z^2;
P=struct('x',[x;y],'f',f);

% Create solver instance
F=nlpsol('F','ipopt',P);

% Solve the problem
r=F('x0',[2.5 3.0])
disp(r.x)
