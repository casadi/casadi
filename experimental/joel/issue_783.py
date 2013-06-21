from casadi import *

x = msym("x",2,1)

fun = MXFunction([x],[diag(x[[1,0]])])
fun.init()

storage2 = []

def flatten(l):
  ret = []
  for i in l:
    ret.extend(i)
  return ret

for f,sym,Function in [(fun,msym,MXFunction),(fun.expand(),ssym,SXFunction)]:
  f.init()
  print Function
  inputss = [sym("i",f.input(i).sparsity()) for i in range(f.getNumInputs()) ]
  aseeds = [[ sym("a",f.output(i).sparsity()) for i in range(f.getNumOutputs()) ]]

  _,_,adjsens = f.eval(inputss,[],aseeds)

  vf = Function(inputss+flatten([aseeds[0]]),flatten([list(adjsens[0])]))
  vf.init()

  # Added to make sure that the same seeds are used for SX and MX
  if Function==MXFunction:
    vf_mx = vf

  inputss2 = [sym("i",vf_mx.input(i).sparsity()) for i in range(vf.getNumInputs()) ]
  aseeds2 = [[ sym("a",vf_mx.output(i).sparsity()) for i in range(vf.getNumOutputs()) ]]
  _,_,adjsens2 = vf.eval(inputss2,[],aseeds2)

  vf2 = Function(inputss2+flatten([aseeds2[0]]),flatten([list(adjsens2[0])]))
  vf2.init()

  offset = 0
  for i in range(vf2.getNumInputs()):
    vf2.setInput(DMatrix(vf2.input(i).sparsity(),range(offset,offset+vf2.input(i).size())),i)
    offset+=vf2.input(i).size()

  vf2.evaluate()

  storage2.append([vf2.getOutput(i) for i in range(vf2.getNumOutputs())])

print "second-order"
for k,(a,b) in enumerate(zip(storage2[0],storage2[1])):
  #print a , " == ", b
  if b.numel()==0 and sparse(a).size()==0: continue
  if a.numel()==0 and sparse(b).size()==0: continue
  if not(sparse(a-b).size()==0):
    raise Exception("At output(%d) : %s <-> %s" % (k,str(a),str(b)))
