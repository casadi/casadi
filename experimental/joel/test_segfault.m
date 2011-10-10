casadi
x = SX("x")
f = SXFunction({{x}},{{x}})
disp 'test "segfault" ok'
