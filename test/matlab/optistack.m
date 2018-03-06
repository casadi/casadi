
import casadi.*

is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;



opti = Opti();




x = opti.variable();
y = opti.variable();
z = opti.variable();

disp(opti.debug.describe(x))

% Test doesn't work if you supply the contents of this file in the terminal,
% as opposed to supplying the name of this file in the terminal
if ~is_octave && ispc
  assert(~isempty(strfind(opti.debug.describe(x),'optistack.m')))
end

opti.minimize((x-1)^2+(y-2)^2+(z-3)^2)
opti.solver('ipopt')
sol = opti.solve();
assert(norm(sol.value(x)-1)<1e-7)
assert(norm(sol.value(y)-2)<1e-7)
assert(norm(sol.value(z)-3)<1e-7)

opti.subject_to()
opti.subject_to({z<=y, y<=x})
sol = opti.solve();

assert(norm(sol.value(x)-2)<1e-7)
assert(norm(sol.value(y)-2)<1e-7)
assert(norm(sol.value(z)-2)<1e-7)

opti.subject_to()
opti.subject_to({z<=y<=x})
sol = opti.solve();

assert(norm(sol.value(x)-2)<1e-7)
assert(norm(sol.value(y)-2)<1e-7)
assert(norm(sol.value(z)-2)<1e-7)

opti.subject_to()
opti.subject_to({z>=y>=x})
sol = opti.solve();

assert(norm(sol.value(x)-1)<1e-4)
assert(norm(sol.value(y)-2)<1e-4)
assert(norm(sol.value(z)-3)<1e-4)

eps = 1e-5;
fun = @(x,y) { 
          y>=2.5;
          y>=1.5;
          2.5>=y;
          1.5>=y;
          y<=2.5;
          y<=1.5;
          2.5<=y;
          1.5<=y;
          y>=x;
          y<=x;
          y<=0;
          y==1.5;
          y==x;
          3<= y <=4;
          0<= y <=1;
          4>= y >=3;
          1>= y >=0;
          3<= x<= y<= 4;
          3<= y<= x<= 4;
          4>= y>= x>= 3;
          4>= x>= y>= 3;};


opti = Opti();

x = opti.variable();
y = opti.variable();
f = (x-1)^2+(y-2)^2;
opti.minimize(f);

opti.solver('ipopt')

data = load('data1.mat');data = data.data;

tests = fun(x, y);


for i = 1:size(tests,1)

  dual_yalmip = data{3,i};
  xs_yalmip = data{1,i};
  ys_yalmip = data{2,i};
  
  opti.subject_to();
  con = tests{i};
  opti.subject_to(con);
  sol = opti.solve();
  assert(norm(sol.value(x)-xs_yalmip)<1e-3);
  assert(norm(sol.value(y)-ys_yalmip)<1e-3);
  dual_1 = sol.value(opti.dual(con));
  assert(norm(dual_1(:)-dual_yalmip)<1e-3);  
end

s = [0.9;1.1;1];

eps = 1e-5;

opti = Opti();

x = opti.variable(3,1);
y = opti.variable(3,1);

f = sum((x-1).^2+((y-2).*s).^2);
opti.minimize(f);
opti.solver('ipopt');

tests = fun(x, y);

data = load('data2.mat');data = data.data;

for i = 1:size(tests,1)
  
  dual_yalmip = data{3,i};
  xs_yalmip = data{1,i};
  ys_yalmip = data{2,i};
  
  opti.subject_to();
  con = tests{i};
  opti.subject_to(con);
  sol = opti.solve();
  assert(norm(sol.value(x)-xs_yalmip)<1e-3);
  assert(norm(sol.value(y)-ys_yalmip)<1e-3);
  dual_1 = sol.value(opti.dual(con));
  assert(norm(dual_1(:)-dual_yalmip(:))<1e-3);  

end

s = [0.9;1.1;1]';

eps = 1e-5;

opti = Opti();

x = opti.variable(1,3);
y = opti.variable(1,3);

opti.solver('ipopt');
f = sum((x-1).^2+((y-2).*s).^2);
opti.minimize(f);

tests = fun(x, y);

data = load('data3.mat');data = data.data;


for i = 1:size(tests,1)
  dual_yalmip = data{3,i};
  xs_yalmip = data{1,i};
  ys_yalmip = data{2,i};

  opti.subject_to();
  con = tests{i};
  opti.subject_to(con);
  sol = opti.solve();
  assert(norm(sol.value(x)-xs_yalmip)<1e-3);
  assert(norm(sol.value(y)-ys_yalmip)<1e-3);
  dual_1 = sol.value(opti.dual(con));
  assert(norm(dual_1(:)-dual_yalmip(:))<1e-3);  

end

s = [0.9 0.7;1.1 0.9;1 1];

eps = 1e-5;

opti = Opti();

x = opti.variable(3,2);
y = opti.variable(3,2);
opti.solver('ipopt')
f = sum(sum((x-1).^2+((y-2).*s).^2));
opti.minimize(f);

fun = @(x,y) { 
          y(:)>=2.5;
          y(:)>=1.5;
          2.5>=y(:);
          1.5>=y(:);
          y(:)<=2.5;
          y(:)<=1.5;
          2.5<=y(:);
          1.5<=y(:);
          y(:)>=x(:);
          y(:)<=x(:);
          y(:)<=0;
          y==1.5;
          y==x;
          3<= y(:) <=4;
          0<= y(:) <=1;
          4>= y(:) >=3;
          1>= y(:) >=0;
          3<= x(:)<= y(:)<= 4;
          3<= y(:)<= x(:)<= 4;
          4>= y(:)>= x(:)>= 3;
          4>= x(:)>= y(:)>= 3;};

tests = fun(x, y);

data = load('data4.mat');data = data.data;

for i = 1:size(tests,1)
  dual_yalmip = data{3,i};
  xs_yalmip = data{1,i};
  ys_yalmip = data{2,i};


  opti.subject_to();
  con = tests{i};
  opti.subject_to(con);
  sol = opti.solve();
  assert(norm(sol.value(x)-xs_yalmip)<1e-3);
  assert(norm(sol.value(y)-ys_yalmip)<1e-3);
  dual_1 = sol.value(opti.dual(con));
  assert(norm(dual_1(:)-dual_yalmip(:))<1e-3);  

end

if ~is_octave
  opti = Opti()

  opti.solver('ipopt')
  x = opti.variable()
  y = opti.variable()

  p = opti.parameter()
  
  opti.minimize(x^2+y^2)

  opti.callback(@(i) evalin('base',['A=' num2str(opti.debug.value(p)) ';']));
  opti.set_value(p, 3);
  sol = opti.solve();
  disp(A)
  assert(A==3);
  opti.set_value(p, 2);
  sol = opti.solve();
  assert(A==2);
  opti.set_value(p, 3);
  sol = opti.solve();
  assert(A==3);

  A = [];
  opti.callback();
  sol = opti.solve();
  assert(isempty(A));

  B = [];

  opti.callback(@(i) evalin('base',['B=' num2str(opti.debug.value(p)) ';']));
  sol = opti.solve();
  assert(B==3);
end


opti = casadi.Opti();
x = opti.variable();
y = opti.variable();
z = opti.variable();
opti.minimize(x^2+y^2+z^2);
opti.subject_to(x+y>=z);
opti.subject_to(-3<=x<=3);
opti.subject_to(-3<=y<=3);
opti.subject_to(x+y>=9);
opti.solver('ipopt');
try
    opti.solve();
catch
    
end

opti.debug.show_infeasibilities;