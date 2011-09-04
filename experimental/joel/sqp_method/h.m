% Inequality constraint function
function rhs = h(x)
global par;

% Simulate
[s,v,m] = rocket(x);
[u,dt] = get_u_and_dt(x);

% Calculate the right hand side
rhs = [par.vmax-v; par.umax-u; u-par.umin];