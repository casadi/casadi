from casadi import *

body = '''
real_t a = sqrt(arg[0]>=0 ? arg[0] : -10*arg[0]);
res[0] = a * arg[1];
res[1] = cos(a);
'''
f = jit('myfcn', 2, 2, body, {'compiler':'shell'})

print f([3, 4])
