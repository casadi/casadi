casadi



function [y]=value(x)
  
end

disp('operators on SX')

x = SX("x")
x**2
sin(x)
acos(x)

disp('operators on SXMatrix')

x = symbolic("x")
x**2
sin(x)
acos(x)

disp('SX typemaps')

x=SX("x")
jacobian(x**5,x)

jacobian({x**5},{x})

jacobian([5],{x})

#jacobian({x**5 x**4},{x})

y=symbolic("y",1,1)
jacobian(y**5,y)

#SXFunction(x,x) deliberate fail
SXFunction({x},{x})
#SXFunction({x;x},{x}) deliberate fail

SXFunction({{x}},{{x}})

SXFunction({{x}},{{x x}})

SXFunction({{x}},{{x;2}})

SXFunction({x},{{x,2}})

SXFunction({x},{{x 2; x x}})


SXFunction({y},{y})

SXFunction({y x},{y x})

disp('MX typemaps')

x=MX("x",1,1)


MXFunction({x},{x x})

MXFunction({x},{x 0})

disp('function usage')
x=SX("x")
f = SXFunction({x},{x^2 sin(x)})
f.init()
f.input(0)
x^2
assert(f.getNumInputs()==1)
assert(f.getNumOutputs()==2)

f.input(0).set([2.3])
f.evaluate()
f.output(0)
f.output(1)
assert(f.output(0)(1).toScalar()==2.3^2)
assert(f.output(1)(1).toScalar()==sin(2.3))

disp(x^2)

f = SXFunction({x},{{x^2 sin(x)}})
f.init()
f.input(0)
assert(f.getNumInputs()==1)
assert(f.getNumOutputs()==1)
f.input(0).set([2.3])
f.evaluate()
f.output(0)

assert(f.output(0)(1).toScalar()==2.3**2)
assert(f.output(0)(2).toScalar()==sin(2.3))

y = symbolic("y",2,2)
f = SXFunction({y},{y**2})
f.init()
f.input(0)
assert(f.getNumInputs()==1)
assert(f.getNumOutputs()==1)
f.input(0).set([1 2; 3 4])
f.evaluate()
f.output(0)
assert(f.output(0)(1,1).toScalar()==1)
assert(f.output(0)(1,2).toScalar()==4)
assert(f.output(0)(2,1).toScalar()==9)
assert(f.output(0)(2,2).toScalar()==16)

disp('dense arrays')
s = DMatrix([ 5 6; 4 9])
full(s)

disp('sparse arrays')

s = DMatrix(3,4,[1,2,1],[0,2,2,3],[0.738,0.1,0.99])

disp('bar')

s.toSparse()

disp('foo')


disp('slicing')
s = DMatrix([ 5 6; 4 9])
assert(s(1,1).toScalar()==5)
assert(s(1,2).toScalar()==6)
assert(s(2,1).toScalar()==4)
assert(s(2,2).toScalar()==9)

q = s([1 2])
assert(q(1).toScalar()==5)
assert(q(2).toScalar()==6)

q = s(1:2,:)
assert(q(1,1).toScalar()==5)
assert(q(1,2).toScalar()==6)
assert(q(2,1).toScalar()==4)
assert(q(2,2).toScalar()==9)

q = s(1:2,2)
assert(q(1,1).toScalar()==6)
assert(q(2,1).toScalar()==9)

% need idx_vector
q = s(1:2,[1 2])
assert(q(1,1).toScalar()==5)
assert(q(1,2).toScalar()==6)
assert(q(2,1).toScalar()==4)
assert(q(2,2).toScalar()==9)

q = s([1 2],1)
assert(q(1).toScalar()==5)
assert(q(2).toScalar()==4)

q = s(1,[1 2])
assert(q(1).toScalar()==5)
assert(q(2).toScalar()==6)

disp('slicing assigment')
s = DMatrix([ 5 6; 4 9])
s(1,1) = 4;
assert(s(1,1).toScalar()==4)

s(:,1) = 78;
assert(s(1,1).toScalar()==78)
assert(s(2,1).toScalar()==78)

disp('MX indexing')

x = MX("x",7,8)
q = x(1:3,:)

x(1,1) = 6;
x(:,1) = 3;

disp('overloaded methods')

X=DMatrix([2 1; 1 4])
inv(X)
result=inv(X)*X;
assert(all(full(result)==full(eye(2))))


disp('Operator overloading')

S =               { {DMatrix([1 2; 3 4]),SX("x"),MX("x",1,1)},
                  {3,symbolic("x",2,2),MX("x",2,2)},
                  {DMatrix(3),symbolic("x",2,2),MX("y",2,2)},
		  {[1 2; 3 4],SX("x"),MX("x",1,1)}
                  };
                  
for i=1:numel(S)
  sc = S{i};
  for j = 1:2
    disp("Here we go")
    s = sc{1}
    z = sc{1+j}
    -s;
    -z;
    z+s;
    s+z;
    s.*z;
    z.*s;
    s-z;
    z-s;
    z./s;
    s./z;
    s.^z;
    z.^s;
  end
end

S = {DMatrix(3),symbolic("x",2,2),SX("x"),MX("x",1,1)};
num = {6,DMatrix(6)}

for i=1:numel(S)
  sc = S{i};
  for j=1:2
    rhs = num{j};
    sc
    rhs
    sc/rhs;
    sc*rhs;
    rhs*sc;
    sc^rhs;
  end
end

x = SX("x");
y = symbolic("y",2,1);

x*y
y*x

num = {DMatrix([1 2; 3 4]),[1 2; 3 4]};
sym = {symbolic("x",2,2),MX("x",2,2)};

for i=1:2
  for j=1:2
    n = num{i};
    s = sym{i};
    n*s;
    s*n;
  end
end

disp("Generic_type")
m=2


is_differential_ivec = IVector(2*m);
is_differential_gentype = GenericType(is_differential_ivec)
assert(is_differential_gentype.isIntVector())

is_differential_ivec = [3,4];
is_differential_gentype = GenericType(is_differential_ivec)
is_differential_gentype.isString()
assert(is_differential_gentype.isDoubleVector())

x=SX("x")
f = SXFunction({x},{x})

integrator = CVodesIntegrator(f)

integrator.setOption('is_differential',[1,3]);


x=symbolic("x",3,4)
size(x)

disp("Issue 145")
t = SX("t")
T = SX("T")
T_dot = SX("T_dot")


ffcn_in = cell(1,DAE_NUM_IN);
ffcn_in{1+DAE_T} = t;
ffcn_in{1+DAE_Y} = T;
ffcn_in{1+DAE_YDOT} = T_dot;

ffcn_in

SXFunction(ffcn_in,{t})


x=symbolic("x",3,4)
size(x)
