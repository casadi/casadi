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
import numpy
import unittest
from types import *
from helpers import *

conics = []
if has_nlpsol("ipopt"):
  ipopt_options = {"fixed_variable_treatment":"relax_bounds",
                   "jac_c_constant":"yes",
                   "jac_d_constant":"yes",
                   "hessian_constant":"yes",
                   "tol":1e-12}
  conics.append(("nlpsol",{"nlpsol":"ipopt", "nlpsol_options.ipopt": ipopt_options},{"quadratic": True}))

if has_nlpsol("worhp"):
  worhp_options = {"TolOpti":1e-13}
  conics.append(("nlpsol",{"nlpsol":"worhp", "nlpsol_options.worhp": worhp_options},{"less_digits":1,"quadratic": True}))


if has_conic("ooqp"):
  conics.append(("ooqp",{},{"less_digits":1,"quadratic": True}))

if has_conic("qpoases"):
  conics.append(("qpoases",{},{"quadratic": True}))

if has_conic("cplex"):
  conics.append(("cplex",{},{"quadratic": True}))

# if has_conic("sqic"):
#   conics.append(("sqic",{},{}))

if has_conic("clp"):
  conics.append(("clp",{"verbose":True},{"quadratic": False}))

if has_conic("activeset"):
  conics.append(("activeset",{},{"quadratic": True}))

print(conics)

class ConicTests(casadiTestCase):

  def test_general_convex_dense(self):
    self.message("Convex dense QP with solvers: " + str([conic for conic,options,aux_options in conics]))
    H = DM([[1,-1],[-1,2]])
    G = DM([-2,-6])
    A =  DM([[1, 1],[-1, 2],[2, 1]])

    LBA = DM([-inf]*3)
    UBA = DM([2, 2, 3])

    LBX = DM([0]*2)
    UBX = DM([inf]*2)

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}

    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      self.message("general_convex: " + str(conic))

      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

      try:
        less_digits=aux_options["less_digits"]
      except:
        less_digits=0

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(**solver_in)
      if conic!="activeset": self.assertTrue(solver.stats()["success"])

      self.assertAlmostEqual(solver_out["x"][0],2.0/3,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],4.0/3,max(1,6-less_digits),str(conic))

      self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      self.checkarray(solver_out["lam_a"],DM([3+1.0/9,4.0/9,0]),str(conic),digits=max(1,6-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-8-2.0/9,max(1,6-less_digits),str(conic))

      solver_in["h"]=H*4

      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],1,max(1,3-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],1,max(1,3-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["cost"],-6,max(1,6-less_digits),str(conic))

      self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      self.checkarray(solver_out["lam_a"],DM([2,0,0]),str(conic),digits=max(1,2-less_digits))

      solver_in["h"]=0

      if 'qcqp' in str(conic): continue # Singular hessian
      if 'activeset' in str(conic): continue # Singular hessian

      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["x"][0],2.0/3,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],4.0/3,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["cost"],-9-1.0/3,max(1,6-less_digits),str(conic))

      self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      self.checkarray(solver_out["lam_a"],DM([10.0/3,4.0/3,0]),str(conic),digits=max(1,4-less_digits))

      solver_in["lba"]=[-inf]*3 #  Upper _and_ lower
      solver_in["uba"]=[inf]*3  #  bounds infinite?

      solver_in["ubx"]=5

      if "worhp" in str(qp_options):
        with self.assertRaises(Exception):
          solver_out = solver(solver_in)
        continue
      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],5,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],5,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["cost"],-40,max(1,5-less_digits),str(conic))

      self.assertAlmostEqual(solver_out["lam_x"][0],2,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["lam_x"][1],6,max(1,6-less_digits),str(conic))

      self.checkarray(solver_out["lam_a"],DM([0,0,0]),str(conic),digits=max(1,4-less_digits))

if __name__ == '__main__':
    unittest.main()
