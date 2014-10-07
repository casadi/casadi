#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
from casadi import *
from casadi.tools import *

class OptimizationObject(MX):
  def create(self,shape,name):
    if not isinstance(shape,tuple): shape = (shape,)
    MX.__init__(self,MX.sym(name,Sparsity.dense(*shape)))
    self.mapping[hash(self)] = self
    
class OptimizationVariable(OptimizationObject):
  mapping = {}
 
  def __init__(self,shape=1,lb=-inf,ub=inf,name="x",init=0):
    self.lb, self.ub, self.init, self.name = lb, ub, init, name
    self.create(shape,name)
    
class OptimizationParameter(OptimizationObject):
  mapping = {}

  def __init__(self,shape=1,value=0,name="p"):
    self.value = value
    self.create(shape,name)

class OptimizationContinousVariable(OptimizationObject):  
  def __init__(self,shape=1,lb=-inf,ub=inf,name="x",init=0):
    self.lb, self.ub, self.init, self.name = lb, ub, init, name
    self.create(shape,name)
    self.start = MX.sym("%s(start)" % name,self.sparsity())
    self.end = MX.sym("%s(end)" % name,self.sparsity())
    self.lim_mapping[hash(self.start)] = self
    self.lim_mapping[hash(self.end)] = self

class OptimizationState(OptimizationContinousVariable):
  mapping = {}
  lim_mapping = {}
  
  def __init__(self,shape=1,lb=-inf,ub=inf,name="x",init=0):
    OptimizationContinousVariable.__init__(self,shape=shape,lb=lb,ub=ub,name=name,init=init)

class OptimizationControl(OptimizationContinousVariable):
  mapping = {}
  lim_mapping = {}
    
  def __init__(self,shape=1,lb=-inf,ub=inf,name="u",init=0):
    OptimizationContinousVariable.__init__(self,shape=shape,lb=lb,ub=ub,name=name,init=init)
    
var, par, state, control = OptimizationVariable, OptimizationParameter, OptimizationState, OptimizationControl

def ocp(f,gl=[],verbose=False,N=20,T=1.0):
  if not isinstance(gl,list): raise Exception("Constraints must be given as a list")
  f = f + OptimizationParameter()
  # Determine nature of constraints, either g(x)<=0 or g(x)==0
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
      
  # Get an exhausive list of all casadi symbols that make up f and gl
  vars = set(getSymbols(veccat([f,T]+gl_pure)))
  
  while True:
    newvars = set() 
    for v in vars:
      if hash(v) in OptimizationState.mapping:
        newvars.update(set(getSymbols(OptimizationState.mapping[hash(v)].dot)))
      for m in [OptimizationState.lim_mapping, OptimizationControl.lim_mapping]:
        if hash(v) in m: newvars.update({m[hash(v)]})
    v0 = len(vars)
    vars.update(newvars)
    if v0==len(vars): break
  
  # Find out which OptimizationParameter and 
  # OptimizationVariable objects correspond to those casadi symbols
  syms = {"w":[],"p":[],"u":[],"x":[]}
  for v in vars:
    for name, cl in zip(["w","p","x","u"],[OptimizationVariable, OptimizationParameter, OptimizationState, OptimizationControl]):
      if hash(v) in cl.mapping: syms[name].append(cl.mapping[hash(v)])

  lims = [i.start for i in syms["x"]] + [i.end for i in syms["x"]] + [i.start for i in syms["u"]]  + [i.end for i in syms["u"]]
  
  # Create structures
  states   = struct_symMX([entry(str(hash(i)),shape=i.sparsity()) for i in syms["x"]])
  controls = struct_symMX([entry(str(hash(i)),shape=i.sparsity()) for i in syms["u"]])
  
  X = struct_symMX([entry(str(hash(i)),shape=i.sparsity()) for i in syms["w"]]+[entry("X",struct=states,repeat=N+1),entry("U",struct=controls,repeat=N)])
  P = struct_symMX([entry(str(hash(i)),shape=i.sparsity()) for i in syms["p"]])
  
  ode_out = MXFunction(syms["x"]+syms["u"]+syms["p"]+syms["w"],[((T+0.0)/N)*vertcat([i.dot for i in syms["x"]])])
  ode_out.init()
  
  nonstates = struct_symMX([entry("controls",struct=controls),entry("p",struct=P)]+[entry(str(hash(i)),shape=i.sparsity()) for i in syms["w"]])
  
  ode = MXFunction(daeIn(x=states,p=nonstates),daeOut(ode=ode_out(states[...]+nonstates["controls",...]+nonstates["p",...]+nonstates[...][2:])[0]))
  ode.init()
  
  intg=explicitRK(ode,1,4)

  h_out = MXFunction(syms["x"]+syms["u"]+syms["p"]+syms["w"],[a for a in gl_pure if dependsOn(a,syms["x"]+syms["u"])])
  g_out = MXFunction(syms["p"]+syms["w"]+lims,[a for a in gl_pure if not dependsOn(a,syms["x"]+syms["u"])])
  f_out = MXFunction(syms["p"]+syms["w"]+lims,[f])
  
  for i in [h_out, g_out, f_out, intg]:i.init()
  
  Pw = P[...]+X[...][:len(syms["w"])]
  Lims = X["X",0,...]+X["X",-1,...]+X["U",0,...]+X["U",-1,...]
  
  # Construct NLP constraints
  G = struct_MX(
    [entry(str(i),expr=g) for i,g in enumerate(g_out(Pw+Lims))] + 
    [entry("path",expr=[ h_out(X["X",k,...]+X["U",k,...]+Pw) for k in range(N)]),
     entry("shooting",expr=[ X["X",k+1] - intg(integratorIn(x0=X["X",k],p=veccat([X["U",k]]+Pw)))[0] for k in range(N)])]
  )
  
  nlp = MXFunction(nlpIn(x=X,p=P),nlpOut(f=f_out(Pw+Lims)[0],g=G))
  nlp.setOption("name","nlp")
  nlp.init()

  # Allocate an ipopt solver
  solver = NlpSolver("ipopt",nlp)
  if not verbose:
    solver.setOption("print_time",False)
    solver.setOption("print_level",0)
    solver.setOption("verbose",False)
  solver.init()

  # Set bounds on variables, set initial value
  x0  = X(solver.input("x0"))
  lbx = X(solver.input("lbx"))
  ubx = X(solver.input("ubx"))

  for i in syms["w"]:
    hs = str(hash(i))
    lbx[hs] = i.lb
    ubx[hs] = i.ub
    x0[hs]  = i.init
    
  for j in "xu":
    for i in syms[j]:
      hs = str(hash(i))
      lbx[j.capitalize(),:,hs] = i.lb
      ubx[j.capitalize(),:,hs] = i.ub
      x0[j.capitalize(),:,hs]  = i.init
    
  # Set parameter values
  par = P(solver.input("p"))
  
  for i in syms["p"]:
    h = str(hash(i))
    par[h] = i.value

  # Set constraint bounds
  lbg = G(solver.input("lbg"))
  ubg = G(solver.input("ubg"))
  
  # Set normal constraints bounds
  for i,eq in enumerate([e for g,e in zip(gl,gl_equality) if not dependsOn(g,syms["x"]+syms["u"])]):
    if eq:
      lbg[str(i)] = ubg[str(i)] = 0
    else:
      lbg[str(i)] = -Inf
      ubg[str(i)] = 0

  # Set path constraints bounds
  for i,eq in enumerate([e for g,e in zip(gl,gl_equality) if dependsOn(g,syms["x"]+syms["u"])]):
    if eq:
      lbg["path",:,i] = ubg["path",:,i] = 0
    else:
      lbg["path",:,i] = -Inf
      ubg["path",:,i] = 0
  
  lbg["shooting",:] = ubg["shooting",:] = 0

  # Solve the problem numerically
  solver.evaluate()
  
  # Raise an exception if not converged
  if solver.getStat('return_status')!="Solve_Succeeded":
    raise Exception("Problem failed to solve. Add verbose=True to see what happened.")
    
  # Add the solution to the OptimizationObjects
  opt = X(solver.output("x"))
  
  # Extract solutions
  for i in syms["w"]: i.sol = opt[str(hash(i))]
  for i in syms["x"]: i.sol = opt["X",:,str(hash(i))]
  for i in syms["u"]: i.sol = opt["U",:,str(hash(i))]    
    
  # Return optimal cost
  return float(solver.getOutput("f"))
