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
    
  def __iter__(self):
    while True:
      yield 123

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
  vars = getSymbols(veccat([f]+gl_pure))
  
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
  G = struct_MX([entry(str(i),expr=g) for i,g in enumerate(gl_pure)])

  # Subsitute the casadi symbols for the structured variants
  original = MXFunction(x+p,nlpOut(f=f,g=G))
  original.init()
  
  nlp = MXFunction(nlpIn(x=X,p=P),original(X[...]+P[...]))
  nlp.init()
  
  # Allocate an ipopt solver
  solver = NlpSolver("ipopt",nlp)
  if not verbose:
    solver.setOption("print_time",False)
    solver.setOption("print_level",0)
    solver.setOption("verbose",False)
  solver.setOption("expand",True)
  solver.init()

  # Set bounds on variables, set initial value
  x0 = X(solver.input("x0"))
  lbx = X(solver.input("lbx"))
  ubx = X(solver.input("ubx"))

  for i in x:
    h = str(hash(i))
    lbx[h] = i.lb
    ubx[h] = i.ub
    x0[h] = i.init
    
  # Set parameter values
  par = P(solver.input("p"))
  
  for i in p:
    h = str(hash(i))
    par[h] = i.value

  # Set constraint bounds
  lbg = G(solver.input("lbg"))
  ubg = G(solver.input("ubg"))
  
  for i,eq in enumerate(gl_equality):
    if eq:
      lbg[str(i)] = ubg[str(i)] = 0
    else:
      lbg[str(i)] = -Inf
      ubg[str(i)] = 0
  
  # Solve the problem numerically
  solver.evaluate()
  
  # Raise an exception if not converged
  if solver.getStat('return_status')!="Solve_Succeeded":
    raise Exception("Problem failed to solve. Add verbose=True to see what happened.")
    
  # Add the solution to the OptimizationObjects
  opt = X(solver.output("x"))
  for i in x:
    i.sol = opt[str(hash(i))]
    
  # Return optimal cost
  return float(solver.getOutput("f"))
  
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
