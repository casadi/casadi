import casadi.*

% SX stuff
x = SX.sym('x');
y = SX.sym('y');
f = cos(x*y)
y = SX.ones(3,2)
f = cos(x*y)


y = SX.sym('y',2,1);
z = MX.sym('z',3);




f = SXFunction({x},{cos(x)})
f.init()

f.setInput(3,0)

f.evaluate()


disp(f.getOutput())

res = f.getOutput()-DMatrix(cos(3))
assert(res.isZero())


x = SX.sym('x',4);

f = SXFunction({x},{x(2),x(IMatrix(2)),x(2,1),x(IMatrix(2),IMatrix(1)),x(2:2),x(2:2,1)});
f.init()

f.setInput([1,2,3,4])
f.evaluate()

for i=1:f.getNumOutputs()
  res = f.getOutput(i-1)-2;

  assert(res.isZero())

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

f = MXFunction({x},{x(2),x(IMatrix(2)),x(2,1),x(IMatrix(2),IMatrix(1)),x(2:2),x(2:2,1)});
f.init()

f.setInput([1,2,3,4])
f.evaluate()

for i=1:f.getNumOutputs()
  res = f.getOutput(i-1)-2;

  assert(res.isZero())

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

clear



