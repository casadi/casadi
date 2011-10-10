casadi
x = SX("x")
f = SXFunction({{x}},{{x x}})
disp 'test "no matching overload" ok'
