import casadi.*

is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

n = 5;
x = MX.sym('x',n,n);
y = MX.sym('y',Sparsity.lower(n));
z = MX.sym('z',n,1);
z2 = MX.sym('z2',n,1);

w = x*3;

w = w*y;

w = sin(w);
w = w./DM(magic(n));

w(1:2) = y(1);

g = {w, norm(w,'fro'), abs(w),z'*z,dot(z,z2),bilin(x,z,z2),rank1(x,0.3,z,z2)};

Xc = {x,y};
for i=1:2
  X = Xc{i};

  M = [vec([X;X])];
  M2 = [X(3,:);X(2,:);X(1:3,:)]+3;
  M3 = vertsplit(X,2);
  M4 = vertsplit(vec(X),2);
  M5 = diagcat(X,2*X);
  g = {g{:} M, M2, M3{1}, M4{1},M5,M5+3,2*X};
end
g = {g{:}, 1./x};


args = symvar(veccat(g{:}));
f_mx = Function('f',args,g);
f_sx = f_mx.expand();

f_mx.export_code('matlab','f_mx_exported.m')
f_sx.export_code('matlab','f_sx_exported.m')
clear f_mx_exported
clear f_sx_exported
rehash


try
  rng(1);
catch
  rand ('state', 1)
end

N = f_mx.n_out;

args_num = {};
for i=1:f_mx.n_in
   a = args{i};
   args_num{i} = sparse(casadi.DM(sparsity(a),rand(nnz(a),1)));
end

f_mx_res = cell(1,1);
f_sx_res = cell(1,1);
f_mx_exported_res = cell(1,1);
f_sx_exported_res = cell(1,1);

[f_mx_res{1:N}] = f_mx(args_num{:});
[f_sx_res{1:N}]  = f_sx(args_num{:});

[f_mx_exported_res{1:N}] = f_mx_exported(args_num{:});
[f_sx_exported_res{1:N}] = f_sx_exported(args_num{:});

for i=1:length(f_mx_res)
    assert(norm(full(f_mx_res{i}-f_sx_res{i}))<1e-12);
    assert(norm(full(f_mx_res{i}-f_mx_exported_res{i}))<1e-12);
    assert(norm(full(f_mx_res{i}-f_sx_exported_res{i}))<1e-12);
end

delete('f_mx_exported.m')
delete('f_sx_exported.m')


linsol = Linsol('solver', 'lapackqr', x.sparsity())

r = linsol.solve(x,y,false)
rt = linsol.solve(x,y,true)

f_mx = Function('f',args,{inv_node(x), inv_node(y), det(x), det(y),x\y,r,rt});
f_mx.disp(true)
f_mx.export_code('matlab','f_mx_exported.m')
clear f_mx_exported
rehash

fref = @(x,y) {inv(x),inv(y),det(x),det(y),x\y,x\y,(x'\y)'};


N = f_mx.n_out;

args={x,y};
args_num = {};
for i=1:numel(args)
   a = args{i};
   args_num{i} = sparse(casadi.DM(sparsity(a),rand(nnz(a),1)));
end

f_mx_res = cell(1,1);
f_mx_exported_res = cell(N,1);

f_mx_res = fref(args_num{:});

[f_mx_exported_res{1:N}] = f_mx_exported(args_num{:});

for i=1:length(f_mx_res)
    assert(norm(full(f_mx_res{i}-f_mx_exported_res{i}))<1e-12);
end

delete('f_mx_exported.m')


if ~is_octave
% Seems like octave does not support this

  A = Sparsity.upper(3);

  save('test.mat','A');

  d = load('test.mat');

  assert(d.A==A);


  A = DM(Sparsity.upper(3),rand(6,1));

  save('test.mat','A');

  d = load('test.mat');

  assert(norm(full(d.A-A))==0);


  x = SX.sym('x',3);
  A = SX.sym('A',Sparsity.lower(3));

  x0 = rand(3,1);
  A0 = sparse(DM(Sparsity.lower(3),rand(6,1)));

  f = Function('f',{x,A},{A*x});

  save('test.mat','f');
  d = load('test.mat');

  f_comp = d.f;

  assert(norm(full(f(x0,A0)-f_comp(x0,A0)))<1e-9);


  if 0
	  x = MX.sym('x',3);
	  A = MX.sym('A',Sparsity.lower(3));

	  x0 = rand(3,1);
	  A0 = sparse(DM(Sparsity.lower(3),rand(6,1)));

	  f = Function('f',{x,A},{A*x});

	  save('test.mat','f');
	  d = load('test.mat');

	  f_comp = d.f;

	  assert(norm(full(f(x0,A0)-f_comp(x0,A0)))<1e-9);
  end

end




