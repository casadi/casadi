% Equality constraint function
function rhs = g(x)
global par;

% Simulate
[s,v,m] = rocket(x);

% Calculate the right hand side
rhs = [s(101)-par.s101; v(101)-par.v101];