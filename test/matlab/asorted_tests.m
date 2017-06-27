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

res = r-DM(cos(3))
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
diary on
x = DM.ones(2,2);

x.print_dense();

x = SX.sym('x');


ode = struct('x', x, 'ode', x)

opts = struct
opts.verbose = true

intg = casadi.integrator('integrator', 'rk', ode, opts);
intg.call(struct);
diary off

logged = fileread('diary');
assert(~isempty(strfind(logged,'1')))


clear

% Are we running Octave?
is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;


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
assert(~isempty(strfind(msg,'  FUNCTION(char,{SX},{SX},struct)')))
assert(~isempty(strfind(msg,'You have: char, SX')))

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


if Importer.has_plugin('clang')
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

if ~is_octave
  [T,X] = ode15s(returntypes('full',ode),[0 300],[2 0]);

  [T,X] = ode15s(returntypes('full',ode),[0 300],[2 0]);

  Jode = Function('ode',{t,x},{jacobian(rhs,x)});
  options = odeset('Jacobian',returntypes('full',Jode));
  [T,X] = ode15s(returntypes('full',ode),[0 300],[2 0],options);
end


% boolvector typemap
x = SX.sym('x');
y = SX.sym('y');
f = (x-2.6)^2 + (y-3.6)^2;
nlp = struct('x', [x;y], 'f', f');

if has_nlpsol('bonmin')
  options = struct;
  options.discrete = [false,false];
  solver = nlpsol('solver', 'bonmin', nlp,options);
  sol = solver('x0',[1 1]);

  assert(all(full(sol.x)==[2.6;3.6]))

  options = struct;
  options.discrete = [true,false];
  solver = nlpsol('solver', 'bonmin', nlp,options);
  sol = solver('x0',[1 1]);

  assert(all(full(sol.x)==[3;3.6]))

  options = struct;
  options.discrete = [false,true];
  solver = nlpsol('solver', 'bonmin', nlp,options);
  sol = solver('x0',[1 1]);

  assert(all(full(sol.x)==[2.6;4]))

  options = struct;
  options.discrete = [true,true];
  solver = nlpsol('solver', 'bonmin', nlp,options);
  sol = solver('x0',[1 1]);

  assert(all(full(sol.x)==[3;4]))
  
  options = struct;
  options.discrete = {true,true};
  solver = nlpsol('solver', 'bonmin', nlp,options);
  sol = solver('x0',[1 1]);

  assert(all(full(sol.x)==[3;4]))
end

data = { [1 3;11 17] , [1 3] , [1;3] 3};

for i=1:numel(data)
  A = data{i};
  B = reshape(DM(A),size(A));
  assert(all(sum(A)==full(sum(B))))
  assert(all(sum(A,1)==full(sum(B,1))))
  assert(all(sum(A,2)==full(sum(B,2))))

  if isvector(A)
    assert(all(norm(A)==full(norm(B))))
    assert(all(norm(A,1)==full(norm(B,1))))
    assert(all(norm(A,2)==full(norm(B,2))))
    assert(all(norm(A,inf)==full(norm(B,inf))))
    assert(all(norm(A,'inf')==full(norm(B,'inf'))))
    assert(all(norm(A,'fro')==full(norm(B,'fro'))))
  else
    assert(all(norm(A,'fro')==full(norm(B,'fro'))))
  end
end

assert(all(size(DM([1 3]))==[1 2]))
assert(all(size(DM([1;3]))==[2 1]))

x=SX.sym('x');
f= Function('f',{x},{x});

f.map('mymap','serial',100,[0],[0])
f.map('mymap','serial',100,{0},{0})

x = SX.sym('x')
y = SX.sym('y',2)
z = SX.sym('z',2,2)
a = SX.sym('a',Sparsity.upper(2))

if ~(is_octave & ismac)
  f = Function('f',{x,y,z,a},{x,y,z,a})
  F = returntypes('full',f);

  [a,b,c,d] = F(1,2,3,4);

  assert(~issparse(c));
  assert(~issparse(d));

  F = returntypes('sparse',f);

  [a,b,c,d] = F(1,2,3,4);

  assert(issparse(c));
  assert(issparse(d));

  F = returntypes('full|sparse',f);

  [a,b,c,d] = F(1,2,3,4);

  assert(~issparse(c));
  assert(issparse(d));

  F = returntypes({'sparse','full','sparse','full'},f);

  [a,b,c,d] = F(1,2,3,4);

  assert(issparse(a));
  assert(~issparse(b));
  assert(issparse(c));
  assert(~issparse(d));
end

Sparsity(3,3)

x = MX.sym('x',3,3);

assert(numel(x)==9);

xr = tril(x);
xnz = xr{:};

assert(numel(xnz)==6);



if ~is_octave 
  x=SX.sym('x');
  y = 4;

  warning('dummy')
  save('test.mat','x','y');
  assert(~isempty(strfind(lastwarn,'not supported')))

  warning('dummy')
  data = load('test.mat');
  assert(~isempty(strfind(lastwarn,'not supported')))
  assert(data.y==4)
  assert(data.x.isnull)
end

Xs = {SX, MX};
for j=1:2;
  X = Xs{j};
  
  for sA=[1,3]
     for sy=[3]
        
        A = X.sym('A',sA,sA);
        y = X.sym('y',sy);
        
        yv = [7;2;4];
        Av = [13 0.2 1;1 9 2;0.1 1 3];
        yv = yv(1:sy);
        Av = Av(1:sA,1:sA);
        F = Function('f',{A,y},{A\y, A\yv, Av\y, A\DM(yv), DM(Av)\y, DM(Av)\DM(yv)});
        out = F.call({Av,yv});
        for i=1:numel(out)
          assert(norm(Av\yv-full(out{i}),1)<=1e-12);
        end
        
        yv = yv';
        y = y';
        F = Function('f',{A,y},{y/A, yv/A, y/Av, DM(yv)/A, y/DM(Av), DM(yv)/DM(Av)});
        out = F.call({Av,yv});
        for i=1:numel(out)
          assert(norm(yv/Av-full(out{i}),1)<=1e-12);
        end
    end
  end


end




x=SX.sym('x');
f=Function('f',{x},{2*x,DM.eye(2)*x});
f.generate('fmex',struct('mex',true));
if is_octave 
mex fmex.c
else
mex fmex.c -largeArrayDims
end
[a,b] = fmex('f',3);
assert(norm(a-6,1)==0);
assert(norm(b-3*eye(2),1)==0);
assert(issparse(a));
assert(issparse(b));
if is_octave 
mex fmex.c -DCASASI_MEX_ALWAYS_DENSE
else
mex fmex.c -largeArrayDims -DCASASI_MEX_ALWAYS_DENSE
end
[a,b] = fmex('f',3);
assert(norm(a-6,1)==0);
assert(norm(b-3*eye(2),1)==0);
assert(~issparse(a));
assert(~issparse(b));
if is_octave 
mex fmex.c -DCASASI_MEX_ALLOW_DENSE
else
mex fmex.c -largeArrayDims -DCASASI_MEX_ALLOW_DENSE
end
[a,b] = fmex('f',3);
assert(norm(a-6,1)==0);
assert(norm(b-3*eye(2),1)==0);
assert(~issparse(a));
assert(issparse(b));





