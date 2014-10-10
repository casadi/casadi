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
    MX.__init__(self,name,sp_dense(*shape))
    self.mapping[hash(self)] = self

class OptimizationContext:
  eval_cache = {}
  
class OptimizationVariable(OptimizationObject):
  mapping = {}
 
  def __init__(self,shape=1,lb=-inf,ub=inf,name="x",init=0):
    """
      Create a decision variable
      
      Parameters
      -------------------
      
      shape: integer or (integer,integer)
        Matrix shape of the symbol
      
      name: string
        A name for the symbol to be used in printing.
        Not required to be unique
   
      lb: number
        Lower bound on the decision variable
        May also be set after initialization as 'x.lb = number'

      ub: number
        Upper bound on the decision variable
        May also be set after initialization as 'x.ub = number'
        
      init: number
        Initial guess for the optimization solver
        May also be set after initialization as 'x.init = number'
        
    """
    self.lb = lb
    self.ub = ub
    self.create(shape,name)
    self.init = init
    self.sol = None
    
class OptimizationParameter(OptimizationObject):
  mapping = {}

  def __init__(self,shape=1,value=0,name="p"):
    """
      Create a parameter, ie a thing that is fixed during optimization
      
      Parameters
      -------------------
      
      shape: integer or (integer,integer)
        Matrix shape of the symbol
      
      name: string
        A name for the symbol to be used in printing.
        Not required to be unique
   
      value: number or matrix
        Value that the parameter should take during optimization
        May also be set after initialization as 'x.value = number'

    """
    self.value = value
    self.create(shape,name)
    
var = OptimizationVariable
par = OptimizationParameter

def maximize(f,**kwargs):
  return - minimize(-f,**kwargs)

def minimize(f,gl=[],verbose=False):
  """
   
   Miminimizes an objective function subject to a list of constraints
   
   Parameters
   -------------------
   
    f:    symbolic expression
       objective function
       
    gl:   list of constraints
       each constraint should have one of these form:
             * lhs<=rhs
             * lhs>=rhs
             * lhs==rhs
             
             where lhs and rhs are expression
             
    Returns
    -------------------
    
    If numerical solution was succesful,
    returns cost at the optimal solution.
    Otherwise raises an exception.
    
    Example
    -------------------
    
    x = var()
    y = var()

    cost = minimize((1-x)**2+100*(y-x**2)**2)
    print "cost = ", cost
    print "sol = ", x.sol, y.sol
  
  """
  
  if not isinstance(gl,list):
    raise Exception("Constraints must be given as a list")
  
  # Determine nature of constraints, either g(x)<=0 or g(x)==0
  g_le = []
  g_eq = []
  g_nsd = []
  for g in gl:
    if g.isOperation(OP_LE) or g.isOperation(OP_LT):
      if (min(g.getDep(0).shape) > 1 and g.getDep(0).shape[0]==g.getDep(0).shape[1]) or (min(g.getDep(1).shape) > 1 and g.getDep(1).shape[0]==g.getDep(1).shape[1]):
        g_nsd.append(g.getDep(0)-g.getDep(1))
      else:
        g_le.append(g.getDep(0)-g.getDep(1))
    elif g.isOperation(OP_EQ):
      g_eq.append(g.getDep(0)-g.getDep(1))
    else:
      print g
      raise Exception("Constrained type unknown. Use ==, >= or <= .")
      
  # Get an exhausive list of all casadi symbols that make up f and gl
  vars = getSymbols(veccat([f]+g_le+g_eq+g_nsd))
  
  
  
  # Find out which OptimizationParameter and 
  # OptimizationVariable objects correspond to those casadi symbols
  x = []
  p = []
  for v in vars:
    if hash(v) in OptimizationVariable.mapping:
      x.append(OptimizationVariable.mapping[hash(v)])
    elif hash(v) in OptimizationParameter.mapping:
      p.append(OptimizationParameter.mapping[hash(v)])
    else:
      raise Exception("Cannot happen")
  
  # Create structures
  X = struct_symMX([entry(str(hash(i)),shape=i.sparsity()) for i in x])
  P = struct_symMX([entry(str(hash(i)),shape=i.sparsity()) for i in p])
  G_le = struct_MX([entry(str(i),expr=g) for i,g in enumerate(g_le)])
  G_eq = struct_MX([entry(str(i),expr=g) for i,g in enumerate(g_eq)])
  G_nsd = struct_MX([entry(str(i),expr=g) for i,g in enumerate(g_nsd)])

  # Subsitute the casadi symbols for the structured variants
  original = MXFunction(x+p,[G_le,G_eq,G_nsd,f])
  original.init()
  
  G = original.evalMX(X[...]+P[...])
  
  def linear(g):
    return jacobian(g,X),substitute(g,X,DMatrix.zeros(X.size))
  
  (A_le,b_le),(A_eq,b_eq),(A_nsd,b_nsd),(A_f,b_f) = map(linear,G)
  
  if A_le.shape[1]==0:
    A_le = DMatrix.zeros(0,X.size)

  if A_eq.shape[1]==0:
    A_eq = DMatrix.zeros(0,X.size)
    
  if A_nsd.shape[1]==0:
    A_nsd = DMatrix.zeros(0,X.size)
    
  f = MXFunction([P],[A_le,b_le,A_eq,b_eq,A_nsd.T,b_nsd,A_f.T,b_f])
  f.init()
  
  A = vertcat([A_le,A_eq])

  if A.shape[1]==0:
    A = DMatrix.zeros(0,X.size)
    
  G_nsd_block = blkdiag(g_nsd)
  
  Fi = []
  for j in range(X.size):
    a = DMatrix(f.output(4)[j,:].sparsity(),1)
    makeDense(a)
    patt = DMatrix(G_nsd_block.sparsity(),a.data())
    makeSparse(patt)
    Fi.append(patt)

  F = vertcat(Fi)
  
  a = DMatrix(f.output(5).sparsity(),1)
  makeDense(a)
  patt = DMatrix(G_nsd_block.sparsity(),a.data())
  makeSparse(patt)
  G = patt
    
  solver = DSdpSolver(sdpStruct(a=A.sparsity(),f=F.sparsity(),g=G.sparsity()))
  if not verbose:
    solver.setOption("_printlevel",0)
  solver.init()
  
  # Set parameter values
  par = P(f.input(0))
  
  for i in p:
    h = str(hash(i))
    par[h] = i.value
  
  f.evaluate()
    
  # Set bounds on variables
  lbx = X(solver.input("lbx"))
  ubx = X(solver.input("ubx"))

  for i in x:
    h = str(hash(i))
    lbx[h] = i.lb
    ubx[h] = i.ub
  

  # Set constraint bounds
  solver.setInput(vertcat([-DMatrix.inf(f.output(1).shape),-f.output(3)]),"lba")
  solver.setInput(vertcat([-f.output(1),-f.output(3)]),"uba")
  
  solver.setInput(vertcat([f.output(0),f.output(2)]),"a")
    
  solver.input("f").set(f.output(4).data())
  solver.input("g").set((-f.output(5)).data())
  
  solver.setInput(dense(f.output(6)),"c")
  
  # Solve the problem numerically
  solver.evaluate()
  
  # Add the solution to the OptimizationObjects
  opt = X(solver.output("x"))
  for i in x:
    i.sol = opt[str(hash(i))]
    
  OptimizationContext.eval_cache = {}
    
  # Return optimal cost
  return float(solver.getOutput("cost"))
  
def value(e,nums={}):
  """
  Evaluates the expression numerically
  
   Parameters
   -------------------
     e: expression to be evaluated
     
     nums: optional dictionary denoting the values of Variables
       if not supplied, the optimal values are assumed
     
  """
  if e in OptimizationContext.eval_cache:
    f,xp = OptimizationContext.eval_cache[e]
  else:
    # Get an exhausive list of all casadi symbols that make up f and gl
    vars = nums.keys() if nums else getSymbols(e)
    
    # Find out which OptimizationParameter and 
    # OptimizationVariable objects correspond to those casadi symbols
    xp = []
    for v in vars:
      if hash(v) in OptimizationVariable.mapping:
        xp.append(OptimizationVariable.mapping[hash(v)])
      elif hash(v) in OptimizationParameter.mapping:
        xp.append(OptimizationParameter.mapping[hash(v)])
      else:
        raise Exception("Cannot happen")
        
    f = MXFunction(xp,[e])
    f.init()
    OptimizationContext.eval_cache[e] = (f,xp)

  for i in range(len(xp)):  
    f.setInput(nums.get(xp[i],xp[i].sol),i)

  f.evaluate()
  
  return f.output()
