from casadi import *
from copy import copy

x = msym("x",2,1)

z=copy(x)
z[[1,0]]*=x

fun = MXFunction([x],[z])
fun.init()

storage2 = []

def flatten(l):
  ret = []
  for i in l:
    ret.extend(i)
  return ret

for f,sym,Function,X in [(fun,msym,MXFunction,MX),(fun.expand(),ssym,SXFunction,SXMatrix)]:
  f.init()
  print Function

  aseeds = sym("a",2)

  _,_,[[adjsens]] = f.eval([X.ones(2)],[],[[aseeds]])

  vf = Function([aseeds],[adjsens])
  vf.init()

  for i in range(vf.getNumInputs()):
    vf.setInput(DMatrix(vf.input(i).sparsity(),range(vf.input(i).size())),i)

  vf.evaluate()

  # Added to make sure that the same seeds are used for SX and MX
  if Function==MXFunction:
    vf_mx = vf

  inputss2 = [sym("i",vf_mx.input(i).sparsity()) for i in range(vf.getNumInputs()) ]
  fseeds2 = [[ sym("f",vf_mx.input(i).sparsity()) for i in range(vf.getNumInputs())]]
  aseeds2 = [[ sym("a",vf_mx.output(i).sparsity()) for i in range(vf.getNumOutputs()) ]]
  res2,fwdsens2,adjsens2 = vf.eval(inputss2,fseeds2,aseeds2)

  vf2 = Function(inputss2+flatten([fseeds2[0]+aseeds2[0]]),list(res2)+flatten([list(fwdsens2[0])+list(adjsens2[0])]))
  vf2.init()

  offset = 0
  for i in range(vf2.getNumInputs()):
    vf2.setInput(DMatrix(vf2.input(i).sparsity(),range(offset,offset+vf2.input(i).size())),i)
    offset+=vf2.input(i).size()

  vf2.evaluate()

  storage2.append([vf2.getOutput(i) for i in range(vf2.getNumOutputs())])


for k,(a,b) in enumerate(zip(storage2[0],storage2[1])):
  #print a , " == ", b
  if b.numel()==0 and sparse(a).size()==0: continue
  if a.numel()==0 and sparse(b).size()==0: continue
  if not(sparse(a-b).size()==0):
    raise Exception("At output(%d) : %s <-> %s" % (k,str(a),str(b)))
