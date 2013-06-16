from casadi import *
n = 12
m = 3

x = msym("x",n,m)
row_offset = range(0,n,n/3)
print "row_offset = ", row_offset

r1,r2,r3 = vertsplit(x,row_offset)

f = MXFunction([x],[r3,r2])
f.init()
f.setInput(range(n*m))

d = (n*m-1)*[0] + [1]
f.setFwdSeed(d)

f.evaluate(1,0)
for i in range(2): print f.output(i)
for i in range(2): print f.fwdSens(i)

print f

gf = f.derivative(0,1)
print gf



