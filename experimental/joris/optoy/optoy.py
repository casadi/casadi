from casadi import *
from casadi.tools import *

class OptimizationVariable(MX):
  mapping = {}
 
  def __init__(self,shape=1,lb=-inf,ub=inf,name="x",init=0):
    self.lb = lb
    self.ub = ub
    if not isinstance(shape,tuple): shape = (shape,) 
    self._sparsity = sp_dense(*shape)
    self.init = init
    self.sol = None
    MX.__init__(self,name,self._sparsity)
    OptimizationVariable.mapping[hash(self)] = self 

class OptimizationParameter(MX):
  mapping = {}

  def __init__(self,shape=1,value=0,name="p"):
    self.value = value
    if not isinstance(shape,tuple): shape = (shape,) 
    self._sparsity = sp_dense(*shape)
    MX.__init__(self,name,self._sparsity)
    OptimizationParameter.mapping[hash(self)] = self
    
var = OptimizationVariable
par = OptimizationParameter

def maximize(f,**kwargs):
  return - minimize(-f,**kwargs)

def minimize(f,gl=[],verbose=False):
  
  # Determine nature of constraints
  gl_pure = []
  gl_equality = []
  for g in gl:
    if g.isOperation(OP_LE) or g.isOperation(OP_LT):
      gl_pure.append(g.getDep(0)-g.getDep(1))
      gl_equality.append(False)
    elif g.isOperation(OP_EQ):
      gl_pure.append(g.getDep(0)-g.getDep(1))
      gl_equality.append(True)
    else:
      raise Exception("Constrained type unknown. Use ==, >= or <= .")
      
  vars = getSymbols(veccat([f]+gl_pure))
  
  x = []
  p = []
  for v in vars:
    if hash(v) in OptimizationVariable.mapping:
      x.append(OptimizationVariable.mapping[hash(v)])
    elif hash(v) in OptimizationParameter.mapping:
      p.append(OptimizationParameter.mapping[hash(v)])
    else:
      raise Exception("Cannot happen")
  
  X = struct_msym([entry(str(hash(i)),shape=i.sparsity()) for i in x])
  P = struct_msym([entry(str(hash(i)),shape=i.sparsity()) for i in p])
  G = struct_MX([entry(str(i),expr=g) for i,g in enumerate(gl_pure)])

  original = MXFunction(x+p,nlpOut(f=f,g=G))
  original.init()
  
  nlp = MXFunction(nlpIn(x=X,p=P),original.evalMX(X[...]+P[...]))
  nlp.init()
  
  solver = IpoptSolver(nlp)
  if not verbose:
    solver.setOption("print_time",False)
    solver.setOption("print_level",0)
    solver.setOption("verbose",False)
  solver.setOption("expand",True)
  solver.init()

  x0 = X(solver.input("x0"))
  lbx = X(solver.input("lbx"))
  ubx = X(solver.input("ubx"))
  
  par = P(solver.input("p"))
  
  for i in x:
    h = str(hash(i))
    lbx[h] = i.lb
    ubx[h] = i.ub
    x0[h] = i.init
  
  for i in p:
    h = str(hash(i))
    par[h] = i.value

  lbg = G(solver.input("lbg"))
  ubg = G(solver.input("ubg"))
  
  for i,eq in enumerate(gl_equality):
    if eq:
      lbg[str(i)] = ubg[str(i)] = 0
    else:
      lbg[str(i)] = -Inf
      ubg[str(i)] = 0
  
  solver.solve()
  
  if solver.getStat('return_status')!="Solve_Succeeded":
    raise Exception("Problem failed to solve. Add verbose=True to see what happened.")
  opt = X(solver.output("x"))
  for i in x:
    i.sol = opt[str(hash(i))]
    
  return solver.output("f")
