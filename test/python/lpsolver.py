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
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

lpsolvers = []
if QpSolver.hasPlugin("ooqp"):
  lpsolvers.append(("qp",{"qp_solver": "ooqp" },False))

if NlpSolver.hasPlugin("ipopt"):
  lpsolvers.append(("qp",{"qp_solver": "nlp", "qp_solver_options":{"nlp_solver":"ipopt"}},False))

if NlpSolver.hasPlugin("ipopt"):
  lpsolvers.append(("qp.nlp.ipopt",{},False))

if QpSolver.hasPlugin("ooqp"):
  lpsolvers.append(("qp.ooqp",{},False))

# try:  
#   QpSolver.loadPlugin("dsdp")
#   def SDPLpSolver(st):
#     return DSdpSolver(sdpStruct(a=st["a"],f=Sparsity.sparse(0,0),g=Sparsity.sparse(0,0)))
#   lpsolvers.append((SDPLpSolver,{},False))
# except:
#  pass

print lpsolvers

class LpSolverTests(casadiTestCase):

  def testboundsviol(self):
    A = DMatrix([[-1,1],[1,1],[1,-2]])
    LBA = DMatrix([ -inf, 2, -inf ])
    UBA = DMatrix([ 1, -inf, 4 ])
    LBX = DMatrix([ -inf, 0 ])
    UBX = DMatrix([ inf, inf ])
    c = DMatrix([ 2.0, 1.0 ])
    
    for lpsolver, lp_options, re_init in lpsolvers:
      self.message("lpsolver: " + str(lpsolver))

      solver = LpSolver(lpsolver,lpStruct(a=A.sparsity()))
      solver.setOption(lp_options)
      solver.init()

      solver.setInput(c,"c")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")
      
    
      with self.assertRaises(Exception):
        solver.evaluate()

    A = DMatrix([[-1,1],[1,1],[1,-2]])
    LBA = DMatrix([ -inf, 2, -inf ])
    UBA = DMatrix([ 1, inf, 4 ])
    LBX = DMatrix([ -inf, 0 ])
    UBX = DMatrix([ inf, -inf ])
    c = DMatrix([ 2.0, 1.0 ])
    
    for lpsolver, lp_options, re_init in lpsolvers:
      self.message("lpsolver: " + str(lpsolver))

      solver = LpSolver(lpsolver,lpStruct(a=A.sparsity()))
      solver.setOption(lp_options)
      solver.init()

      solver.setInput(c,"c")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")
      
    
      with self.assertRaises(Exception):
        solver.evaluate()
        
  def test_bounds(self):
    #  min  2 x +y
    #   x,y
    #
    #  s.t.     bounds on x

    A = DMatrix.sparse(0,2)
    LBX = DMatrix([ -inf, 0 ])
    UBX = DMatrix([ inf, inf ])
    c = DMatrix([ 2.0, 1.0 ])
    
    for lpsolver, lp_options, re_init in lpsolvers:
      self.message("lpsolver: " + str(lpsolver))

      solver = LpSolver(lpsolver,lpStruct(a=A.sparsity()))
      solver.setOption(lp_options)
      solver.init()

      solver.setInput(c,"c")
      solver.setInput([0,0],"lbx")
      solver.setInput([10,10],"ubx")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([0,0]),str(lpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([-2,-1]),str(lpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([]),str(lpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],0,5,str(lpsolver))
      
      if re_init:
        solver = lpsolver(lpStruct(a=A.sparsity()))
        solver.setOption(lp_options)
        solver.init()
     
      solver.setInput(c,"c")
      solver.setInput([-1,3],"lbx")
      solver.setInput([10,10],"ubx")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([-1,3]),str(lpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([-2,-1]),str(lpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([]),str(lpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],1,5,str(lpsolver))
      
      if re_init:
        solver = lpsolver(lpStruct(a=A.sparsity()))
        solver.setOption(lp_options)
        solver.init()
      
      solver.setInput(c,"c")
      solver.setInput([-1,3],"lbx")
      solver.setInput([inf,inf],"ubx")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([-1,3]),str(lpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([-2,-1]),str(lpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([]),str(lpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],1,5,str(lpsolver))
      
      if re_init:
        solver = lpsolver(lpStruct(a=A.sparsity()))
        solver.setOption(lp_options)
        solver.init()
      
      solver.setInput(-c,"c")
      solver.setInput([-10,-10],"lbx")
      solver.setInput([0,0],"ubx")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([0,0]),str(lpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([2,1]),str(lpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([]),str(lpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],0,5,str(lpsolver))

      if re_init:
        solver = lpsolver(lpStruct(a=A.sparsity()))
        solver.setOption(lp_options)
        solver.init()
      
      solver.setInput(-c,"c")
      solver.setInput([-10,-10],"lbx")
      solver.setInput([-1,3],"ubx")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([-1,3]),str(lpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([2,1]),str(lpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([]),str(lpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],-1,5,str(lpsolver))
      
      if re_init:
        solver = lpsolver(lpStruct(a=A.sparsity()))
        solver.setOption(lp_options)
        solver.init()
      
      solver.setInput(-c,"c")
      solver.setInput([-inf,-inf],"lbx")
      solver.setInput([-1,3],"ubx")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([-1,3]),str(lpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([2,1]),str(lpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([]),str(lpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],-1,5,str(lpsolver))
      
  def test_all(self):
    #  min  2 x + y
    #   x,y
    #   
    #   s.t.
    #              -x+y  <= 1
    #         2 <=  x+y
    #               x-2y <= 4  
    #
    #         0 <=     y       
    #

    A = DMatrix([[-1,1],[1,1],[1,-2]])
    LBA = DMatrix([ -inf, 2, -inf ])
    UBA = DMatrix([ 1, inf, 4 ])
    LBX = DMatrix([ -inf, 0 ])
    UBX = DMatrix([ inf, inf ])
    c = DMatrix([ 2.0, 1.0 ])
    
    for lpsolver, lp_options, re_init in lpsolvers:
      self.message("lpsolver: " + str(lpsolver))

      solver = LpSolver(lpsolver,lpStruct(a=A.sparsity()))
      solver.setOption(lp_options)
      solver.init()

      solver.setInput(c,"c")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([0.5,1.5]),str(lpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,0]),str(lpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([0.5,-1.5,0]),str(lpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],2.5,5,str(lpsolver))
      
      if re_init:
        solver = lpsolver(lpStruct(a=A.sparsity()))
        solver.setOption(lp_options)
        solver.init()
      
      # Make a LBX active
      solver.setInput(c,"c")
      solver.setInput(A,"a")
      solver.setInput([-inf,2],"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([1,2]),str(lpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,-3]),str(lpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([2,0,0]),str(lpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],4,5,str(lpsolver))

      if re_init:
        solver = lpsolver(lpStruct(a=A.sparsity()))
        solver.setOption(lp_options)
        solver.init()
        
      # Make a UBX active
      solver.setInput(c,"c")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput([inf,1],"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([1,1]),str(lpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,1]),str(lpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([0,-2,0]),str(lpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],3,5,str(lpsolver))

      if re_init:
        solver = lpsolver(lpStruct(a=A.sparsity()))
        solver.setOption(lp_options)
        solver.init()
      # Make both LBX and UBX active
      solver.setInput(c,"c")
      solver.setInput(A,"a")
      solver.setInput([-inf,3],"lbx")
      solver.setInput([inf,3],"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([2,3]),str(lpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,-3]),str(lpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([2,0,0]),str(lpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],7,5,str(lpsolver))

      if re_init:
        solver = lpsolver(lpStruct(a=A.sparsity()))
        solver.setOption(lp_options)
        solver.init()
      # Linear equality constraint
      solver.setInput(c,"c")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput([1, 2, -inf ],"lba")
      solver.setInput([1,inf,4],"uba")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([0.5,1.5]),str(lpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,0]),str(lpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([0.5,-1.5,0]),str(lpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],2.5,5,str(lpsolver))
      
if __name__ == '__main__':
    unittest.main()
