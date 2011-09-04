function [u,dt]=get_u_and_dt(x)
global par;

switch par.obj
    case 'fuel'
        u = x; dt = par.dt;
    case 'time'
        nx = length(x);
        u = x((1:nx-1)'); dt = x(nx);
end