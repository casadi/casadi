import casadi.*
addpath('defs')

x = MX.sym('x');

c = mycallback();

foo = c.create();

y = foo({x});

f = MXFunction('f',{x},y);

J = f.jacobian();

out = f({5});

out{1}

assert(abs(full(out{1})-25)<1e-12)

out = J({8});
out{1}

assert(abs(full(out{1})-16)<1e-6)

%%

c2 = mycallback2();

foo = c2.create();

y = foo({x});

f = MXFunction('f',{x},y);

J = f.jacobian();

out = f({5});

out{1}

assert(abs(full(out{1})-10)<1e-12)

out = J({8});
out{1}

assert(abs(full(out{1})-2)<1e-6)

%%
c2 = mycallback3();

flag = false;
try
  foo = c2.create();
catch
  flag = true;
end
assert(flag);

% 
c2 = mycallback4(1e-5);

foo = c2.create();

y = foo({x});

f = MXFunction('f',{x},y);

J = f.jacobian();

out = J({5});
out{1}

assert(abs(full(out{1})-10)>=1e-6)
assert(abs(full(out{1})-10)<1e-4)
