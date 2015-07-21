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



%%%%%%%% SVD example

n = 3;
m = 3;
c = mysvd(n,m,true);
foo = c.create();

rng(0);

x = rand(n,m);
X = MX.sym('x',n,m);
Y = foo({X});

f = MXFunction('f',{X},Y);
Js = {};
Js_ = {};
for i=1:3
  Js = {Js{:},f.jacobian(0,i-1)};
  J = Js{i};
  Js_ = {Js_{:},J({x})};
end

[u,s,v] = svd(x);

isequalAbs = @(x,y,tol) ( all(all(abs(full(x)-full(y)) <= tol) ));

for j=Js_
  j = j{1};
  assert(isequalAbs(u,j{2},1e-10));
  assert(isequalAbs(diag(s),j{3},1e-10));
  assert(isequalAbs(v,j{4},1e-10));
end


Js_alts = {};
for w=0:1
  c = mysvd(n,m,false);
  foo = c.create();

  Y = foo({X});

  f = MXFunction('f',{X},Y,struct('ad_weight', w));

  J = f.jacobian(0,1);
  Js = {};
  Js_alt = {};
  for i=1:3
    Js = {Js{:},f.jacobian(0,i-1)};
    J = Js{i};
    Js_alt = {Js_alt{:},J({x})};
  end
  
  Js_alts = {Js_alts{:},Js_alt};

  for j=1:length(Js_)
    j0 = Js_{j};
    j1 = Js_alt{j};
    for i=1:length(j0)
      i0 = j0{i};
      i1 = j1{i};
      assert(isequalAbs(i0,i1,1e-5));
    end
  end
end


for j=1:length(Js_)
  J0 = Js_alts{1};
  J1 = Js_alts{2};
  j0 = J0{j};
  j1 = J1{j};
  for i=1:length(j0)
    i0 = j0{i};
    i1 = j1{i};
    assert(isequalAbs(i0,i1,1e-10));
  end
end 

