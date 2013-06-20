from casadi import *

x = msym("x",2,1)

fun = MXFunction([x],[diag(x[[1,0]])])
fun.init()

storage = []

for f,sym,Function in [(fun,msym,MXFunction),(fun.expand(),ssym,SXFunction)]:
  f.init()
  print Function
  r = DMatrix([0,1])
  q = sym("q",2)

  _,[[fwdsens]],_ = f.eval([r],[[q]],[])

  vf = Function([q],[fwdsens])
  vf.init()
  vf.setInput(r)
  vf.evaluate()
  storage.append([vf.getOutput(i) for i in range(vf.getNumOutputs())])

print "first-order"
for k,(a,b) in enumerate(zip(storage[0],storage[1])):
  #print a , " == ", b
  if b.numel()==0 and sparse(a).size()==0: continue
  if a.numel()==0 and sparse(b).size()==0: continue
  if not(sparse(a-b).size()==0):
    raise Exception("At output(%d) : %s <-> %s" % (k,str(a),str(b)))
