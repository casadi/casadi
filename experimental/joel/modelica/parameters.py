import casadi
ocp = casadi.SymbolicOCP()
#ocp.parseFMI('modelDescription.xml')
ocp.parseFMI('modelDescription.xml',{'sort_equations':False,'eliminate_dependent':False})
ocp.sortType(True) # temporary solution: enables the new sorting
print ocp

x = ocp.variable('x')
x_start = ocp.variable('x_start')
u = ocp.variable('u')
u_cost = ocp.variable('u_cost')

print(x_start.getStart()), " == 1.0"
casadi.updateDependent(ocp)
print(u_cost.getStart()), " == 2.0?"
print(u.getMax()), " == 0.2"
x_start.setStart(2)
print(x_start.getStart()), " == 2.0"
print(x.getStart()), " == 2.0"
print(u.getMax()), " == 0.4"
print(u_cost.getStart()), " == 4.0"
u_cost.setStart(3) # Error since u_cost is a dependent parameter?
x.setStart(4) # Not sure what we want here! What should happen if x_start changes value afterwards?

# JModelica
#from pyjmi import CasadiModel
#from pymodelica import compile_fmux
#jn = compile_fmux("Parameters", "parameters.mop")
#model = CasadiModel(jn)
#import matplotlib.pyplot as plt
#res = model.optimize()
#x = res['x']
#u = res['u']
#time = res['time']
#plt.figure(1)
#plt.clf()
#plt.subplot(2, 1, 1)
#plt.plot(time, x)
#plt.xlabel('$t$')
#plt.ylabel('$x$')
#plt.grid(True)

#plt.subplot(2, 1, 2)
#plt.plot(time, u)
#plt.grid(True)
#plt.xlabel('$t$')
#plt.ylabel('$u$')
#plt.show()
