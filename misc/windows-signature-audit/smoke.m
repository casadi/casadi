% Exercise the released MATLAB MEX and dynamically loaded IPOPT plugin.
import casadi.*
disp(which('casadiMEX'));
x = SX.sym('x');
f = Function('square', {x}, {x*x, jacobian(x*x, x)});
[y, dy] = f(3);
assert(full(y) == 9 && full(dy) == 6);
s = nlpsol('solver', 'ipopt', struct('x', x, 'f', (x-2)^2), ...
    struct('ipopt', struct('print_level', 0), 'print_time', false));
r = s('x0', 0);
stats = s.stats();
assert(stats.success && abs(full(r.x)-2) < 1e-7);
disp('PASS: MATLAB symbolic evaluation, AD, and IPOPT');
