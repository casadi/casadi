from casadi import *

x = msym("x",2,1)
y = msym("y",2,2)

w2=x[:]
w2[1]*=x[0]

#fun = MXFunction([x],[w2,blockcat([[1,MX(1,1)],[x[1],x[0]]])])
fun = MXFunction([x,y],[x])
fun.init()

storage = []
storage2 = []

def flatten(l):
  ret = []
  for i in l:
    ret.extend(i)
  return ret

ndir = 2

for f,sym,Function in [(fun,msym,MXFunction),(fun.expand(),ssym,SXFunction)]:
  f.init()
  print Function
  inputss = [sym("i",f.input(i).sparsity()) for i in range(f.getNumInputs()) ]
  fseeds = [[ sym("f",f.input(i).sparsity()) for i in range(f.getNumInputs())] for d in range(ndir)]
  aseeds = [[ sym("a",f.output(i).sparsity()) for i in range(f.getNumOutputs()) ] for d in range(ndir)]

  res,fwdsens,adjsens = f.eval(inputss,fseeds,aseeds)

  vf = Function(inputss+flatten([fseeds[i]+aseeds[i] for i in range(ndir)]),list(res)+flatten([list(fwdsens[i])+list(adjsens[i]) for i in range(ndir)]))
  vf.init()

  for i in range(vf.getNumInputs()):
    vf.input(i).set(DMatrix(vf.input(i).sparsity(),range(vf.input(i).size())))

  vf.evaluate()

  storage.append([DMatrix(vf.output(i)) for i in range(vf.getNumOutputs())])

  inputss2 = [sym("i",vf.input(i).sparsity()) for i in range(vf.getNumInputs()) ]
  fseeds2 = [[ sym("f",vf.input(i).sparsity()) for i in range(vf.getNumInputs())] for d in range(ndir)]
  aseeds2 = [[ sym("a",vf.output(i).sparsity()) for i in range(vf.getNumOutputs()) ] for d in range(ndir)]
  res2,fwdsens2,adjsens2 = vf.eval(inputss2,fseeds2,aseeds2)

  vf2 = Function(inputss2+flatten([fseeds2[i]+aseeds2[i] for i in range(ndir)]),list(res2)+flatten([list(fwdsens2[i])+list(adjsens2[i]) for i in range(ndir)]))
  vf2.init()

  offset = 0
  for i in range(vf2.getNumInputs()):
    vf2.input(i).set(DMatrix(vf2.input(i).sparsity(),range(offset,offset+vf2.input(i).size())))
    offset+=vf2.input(i).size()

  vf2.evaluate()

  storage2.append([DMatrix(vf2.output(i)) for i in range(vf2.getNumOutputs())])

  #print vf2.output(24)

print "first-order"
for k,(a,b) in enumerate(zip(storage[0],storage[1])):
  #print a , " == ", b
  if b.numel()==0 and sparse(a).size()==0: continue
  if a.numel()==0 and sparse(b).size()==0: continue
  if not(sparse(a-b).size()==0):
    raise Exception("At output(%d) : %s <-> %s" % (k,str(a),str(b)))

print "second-order"
for k,(a,b) in enumerate(zip(storage2[0],storage2[1])):
  #print a , " == ", b
  if b.numel()==0 and sparse(a).size()==0: continue
  if a.numel()==0 and sparse(b).size()==0: continue
  if not(sparse(a-b).size()==0):
    raise Exception("At output(%d) : %s <-> %s" % (k,str(a),str(b)))
