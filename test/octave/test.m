casadi

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

disp('function usage')
x=SX("x")
f = SXFunction({x},{x**2 sin(x)})
f.init()
f.input(0)
assert(f.getNumInputs()==1)
assert(f.getNumOutputs()==2)

f.input(0).set([2.3])
f.evaluate()
f.output(1)
assert(f.output(0)(1)==2.3**2)
assert(f.output(1)(1)==sin(2.3))

f = SXFunction({x},{{x**2 sin(x)}})
f.init()
f.input(0)
assert(f.getNumInputs()==1)
assert(f.getNumOutputs()==1)
f.input(0).set([2.3])
f.evaluate()
f.output(0)

assert(f.output(0)(1)==2.3**2)
assert(f.output(0)(2)==sin(2.3))

y = symbolic("y",2,2)
f = SXFunction({y},{y**2})
f.init()
f.input(0)
assert(f.getNumInputs()==1)
assert(f.getNumOutputs()==1)
f.input(0).set([1 2; 3 4])
f.evaluate()
f.output(0)
assert(f.output(0)(1,1)==1)
assert(f.output(0)(1,2)==4)
assert(f.output(0)(2,1)==9)
assert(f.output(0)(2,2)==16)

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
assert(s(1,1)==5)
assert(s(1,2)==6)
assert(s(2,1)==4)
assert(s(2,2)==9)

q = s([1 2])
assert(q(1)==5)
assert(q(2)==6)

q = s(1:2,:)
assert(q(1,1)==5)
assert(q(1,2)==6)
assert(q(2,1)==4)
assert(q(2,2)==9)

q = s(1:2,2)
assert(q(1,1)==6)
assert(q(2,1)==9)

% need idx_vector
q = s(1:2,[1 2])
assert(q(1,1)==5)
assert(q(1,2)==6)
assert(q(2,1)==4)
assert(q(2,2)==9)

q = s([1 2],1)
assert(q(1)==5)
assert(q(2)==4)

q = s(1,[1 2])
assert(q(1)==5)
assert(q(2)==6)

disp('slicing assigment')
s = DMatrix([ 5 6; 4 9])
s(1,1) = 4;
assert(s(1,1)==4)

s(:,1) = 78;
assert(s(1,1)==78)
assert(s(2,1)==78)

disp('MX indexing')

x = MX("x",7,8)
q = x(1:3,:)

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
    s^z;
    z^s;
  end
end

