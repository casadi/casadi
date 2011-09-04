% Objective function
function rhs = f(x)
global par

[u,dt] = get_u_and_dt(x);
switch par.obj
    case 'fuel'
        rhs = 1/2 * transp(u)*u;
    case 'time'
        rhs = 100*dt;
end
