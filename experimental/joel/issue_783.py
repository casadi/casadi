from casadi import *

x = msym("x",2,1)

fun = MXFunction([x],[diag(x[[1,0]])])
fun.init()

storage = []

for f,sym,Function in [(fun,msym,MXFunction),(fun.expand(),ssym,SXFunction)]:
  f.init()
  print Function
  inputss = [sym("i",f.input(i).sparsity()) for i in range(f.getNumInputs()) ]
  fseeds = [[ sym("f",f.input(i).sparsity()) for i in range(f.getNumInputs())]]
  aseeds = [[ sym("a",f.output(i).sparsity()) for i in range(f.getNumOutputs()) ]]

  res,fwdsens,adjsens = f.eval(inputss,fseeds,aseeds)

  vf = Function(inputss+fseeds[0]+aseeds[0],list(res)+list(fwdsens[0])+list(adjsens[0]))
  vf.init()

  for i in range(vf.getNumInputs()):
    vf.setInput(DMatrix(vf.input(i).sparsity(),range(vf.input(i).size())),i)

  vf.evaluate()

  # Added to make sure that the same seeds are used for SX and MX
  if Function==MXFunction:
    vf_mx = vf

  storage.append([vf.getOutput(i) for i in range(vf.getNumOutputs())])

print "first-order"
for k,(a,b) in enumerate(zip(storage[0],storage[1])):
  #print a , " == ", b
  if b.numel()==0 and sparse(a).size()==0: continue
  if a.numel()==0 and sparse(b).size()==0: continue
  if not(sparse(a-b).size()==0):
    raise Exception("At output(%d) : %s <-> %s" % (k,str(a),str(b)))
