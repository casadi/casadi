import casadi.*

x = SX.sym('x');

f = Function('f',{x},{2*x});

F = f.map(2,'thread')

out = F([3 5]);

assert(norm(full(out)-[6 10])==0);
