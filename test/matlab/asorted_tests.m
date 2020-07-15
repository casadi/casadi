import casadi.*

MX(2,1)

% SX stuff
x = SX.sym('x');
y = SX.sym('y');
f = cos(x*y)
y = SX.ones(3,2)
f = cos(x*y)


y = SX.sym('y',2,1);
z = MX.sym('z',3);

DM(true)
SX(true)

r = which_depends([x;y],[x;y])
assert(islogical(r))

f = Function('f',{}, {0})
f = Function('f',cell(1,0), {0})
f = Function('f',cell(0,1), {0})

f = Function('f',{x},{cos(x)})
r = f(3)

disp(r)

res = r-DM(cos(3))
assert(is_zero(res))


x = SX.sym('x',4);

f = Function('f',{x},{x(2),x(DM(2)),x(2,1),x(DM(2),DM(1)),x(2:2),x(2:2,1)});
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

f = Function('f',{x},{x(2),x(DM(2)),x(2,1),x(DM(2),DM(1)),x(2:2),x(2:2,1)});
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
nlp = Function('nlp', {x,p}, {f,g}, char('x','p'), char('f','g'));

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

data = { [1 3 2;11 17 4], [1 3;11 17] , [1 3] , [1;3] 3};

for i=1:numel(data)
  A = data{i};
  B_DM = reshape(DM(A),size(A));
  B_MX = MX(B_DM);
  Bs = {B_DM, B_MX};
  for j=1:2
    B = Bs{j};
    assert(all(sum(A)==full(evalf(sum(B)))))
    assert(all(sum(A,1)==full(evalf(sum(B,1)))))
    assert(all(sum(A,2)==full(evalf(sum(B,2)))))

    assert(all(all(cumsum(A)==full(evalf(cumsum(B))))))
    assert(all(all(cumsum(A,1)==full(evalf(cumsum(B,1))))))
    assert(all(all(cumsum(A,2)==full(evalf(cumsum(B,2))))))

    assert(all(diff(A)==full(evalf(diff(B)))))
    assert(all(diff(A,1)==full(evalf(diff(B,1)))))
    if ~is_octave
      assert(all(diff(A,2)==full(evalf(diff(B,2)))))
    end
    assert(all(all(diff(A,1,1)==full(evalf(diff(B,1,1))))))
    assert(all(all(diff(A,1,2)==full(evalf(diff(B,1,2))))))
    assert(all(all(diff(A,2,1)==full(evalf(diff(B,2,1))))))
    assert(all(all(diff(A,2,2)==full(evalf(diff(B,2,2))))))

    if j==1
      if isvector(A)
        assert(all(norm(A)==full(evalf(norm(B)))))
        assert(all(norm(A,1)==full(evalf(norm(B,1)))))
        assert(all(norm(A,2)==full(evalf(norm(B,2)))))
        assert(all(norm(A,inf)==full(evalf(norm(B,inf)))))
        assert(all(norm(A,'inf')==full(evalf(norm(B,'inf')))))
        assert(all(norm(A,'fro')==full(evalf(norm(B,'fro')))))
      else
        assert(all(norm(A,'fro')==full(evalf(norm(B,'fro')))))
      end
    end
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

%if ~(is_octave & ismac)
if false
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

if isunix
  compile_flags = 'CFLAGS=$CFLAGS -pedantic -std=c89 -fPIC -shared -Wall -Werror -Wextra -Wno-unknown-pragmas -Wno-long-long -Wno-unused-parameter'
else
  compile_flags = '';
end

x=SX.sym('x');
f=Function('f',{x},{2*x,DM.eye(2)*x});
f.generate('fmex',struct('mex',true));
clear fmex
if is_octave
mex -DMATLAB_MEX_FILE fmex.c
else
mex('-largeArrayDims',compile_flags,'fmex.c')
end
[a,b] = fmex('f',3);
assert(norm(a-6,1)==0);
assert(norm(b-3*eye(2),1)==0);
assert(~issparse(a));
assert(issparse(b));
clear fmex
if is_octave
mex -DCASADI_MEX_NO_SPARSE -DMATLAB_MEX_FILE fmex.c
else
mex('-DCASADI_MEX_NO_SPARSE',compile_flags,'-largeArrayDims','fmex.c')
end
[a,b] = fmex('f',3);
assert(norm(a-6,1)==0);
assert(norm(b-3*eye(2),1)==0);
assert(~issparse(a));
assert(~issparse(b));


f.generate('fmex_rp',struct('mex',true,'casadi_int','int','casadi_real','float'));
clear fmex_rp
if is_octave
mex -DMATLAB_MEX_FILE fmex_rp.c
else
mex -largeArrayDims fmex_rp.c
end
[a,b] = fmex_rp('f',3);
assert(norm(a-6,1)==0);
assert(norm(b-3*eye(2),1)==0);
assert(~issparse(a));
assert(issparse(b));
clear fmex_rp
if is_octave
mex -DCASADI_MEX_NO_SPARSE -DMATLAB_MEX_FILE fmex_rp.c
else
mex('-DCASADI_MEX_NO_SPARSE',compile_flags,'-largeArrayDims','fmex_rp.c')
end
[a,b] = fmex_rp('f',3);
assert(norm(a-6,1)==0);
assert(norm(b-3*eye(2),1)==0);
assert(~issparse(a));
assert(~issparse(b));


x=SX.sym('x');
f=Function('f',{x},{2*x,DM.eye(2)*x});
f.generate('fmex',struct('mex',true));
clear fmex
if is_octave
mex -DMATLAB_MEX_FILE fmex.c
else
mex('-largeArrayDims',compile_flags,'fmex.c')
end
[a,b] = fmex(3);
assert(norm(a-6,1)==0);
assert(norm(b-3*eye(2),1)==0);

x=SX.sym('x');
f=Function('f',{x},{2*x,DM.eye(2)*x});
g=Function('g',{x},{3*x,DM.eye(2)*x});
cg = CodeGenerator('fmex', struct('mex',true))
cg.add(f)
cg.add(g)
cg.generate()

clear fmex
if is_octave
mex -DMATLAB_MEX_FILE fmex.c
else
mex('-largeArrayDims',compile_flags,'fmex.c')
end

flag = false;
try
  [a,b] = fmex(3);
catch
  flag = true;
end
assert(flag);

[a,b] = fmex('f',3);
assert(norm(a-6,1)==0);
assert(norm(b-3*eye(2),1)==0);


flag = false;
try
  [a,b] = fmex('h',3);
catch
  flag = true;
end
assert(flag);

flag = false;
try
  [a,b] = fmex('hhhhhh',3);
catch
  flag = true;
end
assert(flag);

x = SX.sym('x',2,2);
y = SX.sym('y',Sparsity.lower(2))
f=Function('f',{x,y},{x+y,3*y});
f.generate('fmex',struct('mex',true));
clear fmex
if is_octave
mex -DMATLAB_MEX_FILE fmex.c
else
mex('-largeArrayDims',compile_flags,'fmex.c')
end
xnom = [1 2;3 5]
ynom = sparse([7 0;6 9])
[a,b] = fmex('f',xnom,ynom);
assert(norm(a-(xnom+ynom),1)==0);
assert(norm(b-3*ynom,1)==0);
assert(~issparse(a));
assert(issparse(b));
xnom = [4 4;4 4]
ynom = sparse([5 0;5 5])
[a,b] = fmex('f',4,5);
assert(norm(a-(xnom+ynom),1)==0);
assert(norm(b-3*ynom,1)==0);
assert(~issparse(a));
assert(issparse(b));

ynom = sparse([5 0;0 6])
[a,b] = fmex('f',4,ynom);
assert(norm(a-(xnom+ynom),1)==0);
assert(norm(b-3*ynom,1)==0);
assert(~issparse(a));
assert(issparse(b));

ynom = sparse([5 0;6 7])
[a,b] = fmex('f',4,sparse([5 12;6 7]));
assert(norm(a-(xnom+ynom),1)==0);
assert(norm(b-3*ynom,1)==0);
assert(~issparse(a));
assert(issparse(b));

ynom = sparse([5 0;6 7])
[a,b] = fmex('f',4,[5 12;6 7]);
assert(norm(a-(xnom+ynom),1)==0);
assert(norm(b-3*ynom,1)==0);
assert(~issparse(a));
assert(issparse(b));

xnom = [4 4;4 4]
ynom = sparse([5 0;5 5])

flag = false;
try
  [a,b] = fmex('f',[2 3],ynom);
catch
  flag = true;
end
assert(flag);

flag = false;
try
  [a,b] = fmex('f',xnom,sparse(ones(3,2)));
catch
  flag = true;
end
assert(flag);

clear fmex
if is_octave
mex -DCASADI_MEX_NO_SPARSE -DMATLAB_MEX_FILE fmex.c
else
mex('-DCASADI_MEX_NO_SPARSE',compile_flags,'-largeArrayDims','fmex.c')
end
xnom = [1 2;3 5]
ynom = sparse([7 0;6 9])
[a,b] = fmex('f',xnom,ynom);
assert(norm(a-(xnom+ynom),1)==0);
assert(norm(b-3*ynom,1)==0);
assert(~issparse(a));
assert(~issparse(b));
xnom = [4 4;4 4]
ynom = sparse([5 0;5 5])
[a,b] = fmex('f',4,5);
assert(norm(a-(xnom+ynom),1)==0);
assert(norm(b-3*ynom,1)==0);
assert(~issparse(a));
assert(~issparse(b));


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

  A = X.sym('A',3,3);
  Av = [13 0.2 1;1 9 2;0.1 1 3];
  for N=[-4,-3,-2,-1,0,1,2,3,4]
    F = Function('f',{A},{A^N,DM(Av)^N});
    out = F.call({Av});
    for i=1:numel(out)
      assert(norm(Av^N-full(out{i}),1)<=1e-12);
    end
  end

end

c = {1,2,4;3 5 6}
e = full(casadi.blockcat(c))
assert(norm(cell2mat(c)-e)==0)

e = full(casadi.blockcat({{1,2,4};{3 5 6}}))
assert(norm(cell2mat(c)-e)==0)
e = full(casadi.blockcat({{1,2,4} {3 5 6}}))
assert(norm(cell2mat(c)-e)==0)
e = full(casadi.blockcat({{1;2;4} {3 5 6}}))
assert(norm(cell2mat(c)-e)==0)
e = full(casadi.blockcat({{1;2;4} {3;5;6}}))
assert(norm(cell2mat(c)-e)==0)

c = casadi.blockcat({1 2 3})
assert(norm(size(c)-[1 3])==0)
c = casadi.blockcat({1;2;3})
assert(norm(size(c)-[3 1])==0)


if ~is_octave
  msg = help('MX/repsum');
  assert(~isempty(strfind(msg,'Given a repeated matrix, computes the sum of repeated parts.')))
end

x = SX.sym('x');
s = StringSerializer();
s.pack(x);
s.pack(sin(x));
 
data = s.encode();
 
s = StringDeserializer(data);
a = s.unpack();
b = s.unpack();

a = sparse([1 0 2; 3 0 4]);
A = DM(a);
assert(nnz(a(:))==nnz(A(:)));
assert(full(norm(a(:)-A(:)))==0);

a(:) = 3;
A(:) = 3;
assert(nnz(a(:))==nnz(A(:)));
assert(full(norm(a(:)-A(:)))==0);

a = sparse([1 0 2; 3 0 4]);
A = DM(a);

r = [1 2 6 7 8 9];
a(:) = r;
A(:) = r;

assert(nnz(a(:))==nnz(A(:)));
assert(full(norm(a(:)-A(:)))==0);

a = sparse([1 0 2; 3 0 4]);
A = DM(a);

r = sparse([0 2 6 0 8 9]);
a(:) = r;
A(:) = r;

assert(nnz(a(:))==nnz(A(:)));
assert(full(norm(a(:)-A(:)))==0);


assert(size(MX(ones(1,4)),1)==1)
assert(size(MX(ones(1,4)),2)==4)



x = MX.sym('x');
y = MX.sym('y');

f = Function('f',{x,y},{x+y},struct('is_diff_in',[1 0]))
f = Function('f',{x,y},{x+y},struct('is_diff_in',[true false]))
f = Function('f',{x,y},{x+y},struct('is_diff_in',{{true false}}))
f = Function('f',{x,y},{x+y},struct('is_diff_in',{{1 0}}))
f = Function('f',{x,y},{x+y},struct('is_diff_in',[1;0]))
f = Function('f',{x,y},{x+y},struct('is_diff_in',[true;false]))
f = Function('f',{x,y},{x+y},struct('is_diff_in',{{true;false}}))
f = Function('f',{x,y},{x+y},struct('is_diff_in',{{1;0}}))


flag = false;
try
  f = Function('f',{x,y},{x+y},struct('is_diff_in',[1 3]))
catch
  flag = true;
end
assert(flag);

f.map(3,[1 0],[0])
f.map(3,[true false],[0])
f.map(3,{true false},[0])
%f.map(3,{1 0},[0])
f.map(3,[1;0],[0])
f.map(3,[true;false],[0])
f.map(3,{true;false},[0])
%f.map(3,{1;0},[0])

flag = false;
try
  f.map(3,[1 3])
catch
  flag = true;
end
assert(flag);


A0 = DM([1 2; 3 4; 5 6]);
A = MX.sym('A',3,2);
x = MX.sym('x');
xi = 1;

%f = Function('f',{A,x},{A{x}});assert(f(A0,xi)==A0{xi})
f = Function('f',{A,x},{A(x,:)});assert(full(norm(f(A0,xi)-A0(xi,:)))==0)
f = Function('f',{A,x},{A(:,x)});assert(full(norm(f(A0,xi)-A0(:,xi)))==0)
f = Function('f',{A,x},{A(1:2,x)});assert(full(norm(f(A0,xi)-A0(1:2,xi)))==0)
f = Function('f',{A,x},{A(x,1:2)});assert(full(norm(f(A0,xi)-A0(xi,1:2)))==0)
f = Function('f',{A,x},{A(x,x)});assert(full(norm(f(A0,xi)-A0(xi,xi)))==0)
f = Function('f',{A,x},{A(x)});assert(full(norm(f(A0,xi)-A0(xi)))==0)

A0 = DM([1 2 3]);
A = MX.sym('A',3);
f = Function('f',{A,x},{A(x)});assert(full(norm(f(A0,xi)-A0(xi)))==0)

A0 = DM([1 2 3])';
A = MX.sym('A',1,3);
f = Function('f',{A,x},{A(x)});assert(full(norm(f(A0,xi)-A0(xi)))==0)


