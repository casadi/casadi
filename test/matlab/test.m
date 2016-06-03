import casadi.*

% SX stuff
x = SX.sym('x');
y = SX.sym('y');
f = cos(x*y)
y = SX.ones(3,2)
f = cos(x*y)


y = SX.sym('y',2,1);
z = MX.sym('z',3);




f = Function('f',{x},{cos(x)})
r = f(3)

disp(r)

res = r{1}-DM(cos(3))
assert(is_zero(res))


x = SX.sym('x',4);

f = Function('f',{x},{x(2),x(IM(2)),x(2,1),x(IM(2),IM(1)),x(2:2),x(2:2,1)});
r = f.call({[1,2,3,4]});

for i=1:f.n_out()
  res = r{i}-2;

  assert(res.is_zero())

end

flag = false;
try
  x(5);
catch
  flag = true;
end
assert(flag);



flag = false;
try
  x(0);
catch
  flag = true;
end
assert(flag);


x = MX.sym('x',4);

f = Function('f',{x},{x(2),x(IM(2)),x(2,1),x(IM(2),IM(1)),x(2:2),x(2:2,1)});
r = f.call({[1,2,3,4]});

for i=1:f.n_out()
  res = r{i}-2;

  assert(res.is_zero())

end


flag = false;
try
  x(5);
catch
  flag = true;
end
assert(flag);

flag = false;
try
  x(0);
catch
  flag = true;
end
assert(flag);



% Checking stdout

delete 'diary'
diary ON
x = DM.ones(2,2);

x.print_dense();

x = SX.sym('x');


ode = struct('x', x, 'ode', x)

opts = struct
opts.verbose = true

intg = casadi.integrator('integrator', 'rk', ode, opts);
intg.call(struct);
diary OFF

logged = fileread('diary');
assert(~isempty(strfind(logged,'1')))


clear


x = SX.sym('x',4,5);
dim = size(x);
assert(dim(1)==4)
assert(dim(2)==5)

% error message beautification

msg = '';
try
  Function('f', [SX.sym('12')])
catch err
  msg = err.message;
end

% See issue #1483
%assert(~isempty(strfind(msg,'  Function(char,{SX} ,{SX} ,Dict)')))
%assert(~isempty(strfind(msg,'You have: char, SX')))

% Check mixing DM and MX
res = (DM(1)+MX(1)) - (MX(1)+DM(1))
assert(is_zero(res))

% Try substitute (non-member function)
a = SX.sym('a');
b = substitute({sin(a), cos(a)}, {a},{a+3});
b{1}
b{2}

veccat(MX.sym('x',2,3),MX.sym('x',4,5))
length(MX.sym('x',4,1))
assert(length(MX.sym('x',4,1))==4)
assert(length(MX.sym('x',1,3))==3)
assert(length(MX.sym('x',4,5))==5)

[n,m] = size(MX.sym('x',4,5));
assert(n==4)
assert(m==5)

% Create function with custom IO scheme
x = SX.sym('x');
p = SX.sym('p');
f = x^2;
g = log(x)-p;
opts = struct('input_scheme', char('x','p'),...
              'output_scheme', char('f','g'));
nlp = Function('nlp', {x,p}, {f,g}, opts);

% Evaluate with numbered inputs and outputs
[res_vec{1:2}] = nlp(1.1, 3.3);

% Evaluate with named inputs and outputs
res_struct = nlp('x',1.1,'p',3.3);
assert(is_zero(res_vec{1}-res_struct.f))
assert(is_zero(res_vec{2}-res_struct.g))

u = SX.sym('u',1,5);

assert(all(size(u(1:2))==[1 2]));
assert(all(size(u([1:2]'))==[1 2]));

u = SX.sym('u',5,1);

assert(all(size(u(1:2))==[2 1]));
assert(all(size(u([1:2]'))==[2 1]));

u = SX.sym('u',4,5);

assert(all(size(u(:,2))==[4 1]));
assert(all(size(u(:,1:2))==[4 2]));

assert(all(size(u(2,:))==[1 5]));
assert(all(size(u(1:2,:))==[2 5]));

assert(all(size(u(1:2,1:3))==[2 3]));


if Importer.hasPlugin('clang')
  x = MX.sym('x');
  F = Function('f',{x},{x^2},struct('jit',true));

  out = F(5);
  assert(full(out)==25)
end


a = DM.ones(Sparsity.upper(5));
i = full(full(sparse(a))-a);
assert(any(any(i))==0);

% Element assignment
A = DM.ones(1,4);
A(2) = 20;
A(3:4) = [30,40];
assert(all(full(A)==[1 20 30 40]))

% Full/sparse Function

x = SX.sym('x',2);
t = SX.sym('t');

rhs = [x(2);1000*(1 - x(1)^2)*x(2) - x(1)];
ode = Function('ode',{t,x},{rhs});

[T,X] = ode15s(full(ode),[0 300],[2 0]);

[T,X] = ode15s(sparse(ode),[0 300],[2 0]);

Jode = Function('ode',{t,x},{jacobian(rhs,x)});
options = odeset('Jacobian',full(Jode));
[T,X] = ode15s(full(ode),[0 300],[2 0],options);
