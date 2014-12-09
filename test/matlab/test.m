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
clear
