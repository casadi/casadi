casadi

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

123

SXFunction({{x}},{{x x}})

SXFunction({{x}},{{x;2}})

SXFunction({x},{{x,2}})

SXFunction({x},{{x 2; x x}})


SXFunction({y},{y})


SXFunction({y x},{y;x})
