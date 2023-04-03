#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
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

class Feasiblesqpmethod:#(casadiTestCase):

  def test_feas_iter_test(self):
      # Set options
      opts = {}
      opts["expand"] = True
      opts["hessian_approximation"] = "exact"
      opts["max_inner_iter"] = 30
      opts["max_iter"] = 20
      opts["feas_tol"] = 1e-6
      opts["tr_scale_vector"] = [1.0, 1.0] #[2.5,3.0,0.75}
      # opts["derivative_test"] = "second-order"

      # Specify QP solver
      opts["qpsol"]  = "nlpsol"
      opts["qpsol_options.nlpsol"] = "ipopt"
      opts["qpsol_options.error_on_fail"] = False
      opts["qpsol_options.nlpsol_options.ipopt.print_level"] = 0
      opts["qpsol_options.nlpsol_options.ipopt.sb"] = "yes"
      opts["qpsol_options.nlpsol_options.print_time"] = 0

      x = SX.sym("x")
      y = SX.sym("y")
      f = x**2 + (y-2)**2
      g = x**2 + y**2
      nlp = {"x": vertcat(x,y), "f": f, "g": g}

      # Create an NLP solver
      # Function solver = nlpsol("solver", "sqpmethod", nlp, opts)
      solver = nlpsol("solver", "feasiblesqpmethod", nlp, opts)

      # Solve the problem
      arg = {}
      arg["x0"] = [2.0, 0.0] #[2.5,3.0,0.75}
      arg["lbg"] = 4
      arg["ubg"] = 4
      res = solver(**arg)
      self.assertEqual(solver.stats()["iter_count"],9)
      self.checkarray(res["x"],vertcat(0,2),digits=7)

  def test_local_conv_test(self):
      """
      The problem contains a problem that contains a region where the
      local SQP method converges and a region from where it diverges. I.e., a 
      globalization strategy is necessary.
      """
      # Set options
      opts = {}
      opts["expand"] = True
      # opts["max_iter"] = 10)
      # opts["verbose"] = true
      # opts["linear_solver"] = "ma57"
      opts["hessian_approximation"] = "exact"
      # opts["derivative_test"] = "second-order"

      # Specify QP solver
      opts["qpsol"]  = "nlpsol"
      opts["qpsol_options.nlpsol"] = "ipopt"
      opts["qpsol_options.error_on_fail"] = False
      opts["qpsol_options.nlpsol_options.ipopt.print_level"] = 0
      opts["qpsol_options.nlpsol_options.ipopt.sb"] = "yes"
      opts["qpsol_options.nlpsol_options.print_time"] = 0

      x = SX.sym("x")
      f = log(exp(x) + exp(-x))
      g = x
      nlp = {"x": x, "f": f, "g": g}

      # Create an NLP solver
      # Function solver = nlpsol("solver", "sqpmethod", nlp, opts)
      solver = nlpsol("solver", "feasiblesqpmethod", nlp, opts)

      print("Here goes the non-converging initial point for local SQP")
      # Solve the problem
      arg = {}
      arg["x0"] = 1.4
      arg["lbg"] = -2
      arg["ubg"] = 2
      res = solver(**arg)
      self.assertEqual(solver.stats()["iter_count"],4)
      self.checkarray(res["x"],[0],digits=3)
      #uout() << "Optimal cost:                     " << double(res.at("f")) << std.endl
      #uout() << "Primal solution:                  " << std.vector<double>(res.at("x")) << std.endl
      #uout() << "Here goes the converging initial point for local SQP" << std.endl
      arg["x0"] = 1.0 #[2.5,3.0,0.75}
      res = solver(**arg)
      self.checkarray(res["x"],[0],digits=4)
      self.assertEqual(solver.stats()["iter_count"],2)

      #  Print solution
      #uout() << "Optimal cost:                     " << double(res.at("f")) << std.endl
      #uout() << "Primal solution:                  " << std.vector<double>(res.at("x")) << std.endl

  def test_qp_test(self):
      """
      The problem is a simple QP. Any SQP solver should just need one
      iteration! Depending on the size of the trust-region. It could take more 
      iteration. But if the trust-region is big enough such that the solution is
      included. The problem should just take one iteration!
      """
      for tr_rad0, iter_count in [(1.0,4),(100.0,2)]:
          # Set options
          opts = {}
          opts["expand"] = True
          # opts["max_iter"] = 10)
          # opts["verbose"] = true
          # opts["linear_solver"] = "ma57"
          opts["hessian_approximation"] = "exact"
          opts["tr_rad0"] = tr_rad0
          # opts["derivative_test"] = "second-order"

          # Specify QP solver
          opts["qpsol"]  = "nlpsol"
          opts["qpsol_options.nlpsol"] = "ipopt"
          opts["qpsol_options.error_on_fail"] = False
          opts["qpsol_options.nlpsol_options.ipopt.print_level"] = 0
          opts["qpsol_options.nlpsol_options.ipopt.sb"] = "yes"
          opts["qpsol_options.nlpsol_options.print_time"] = 0

          x = SX.sym("x")
          f = 2.0*(x+4)**2+2
          g = x
          nlp = {"x": x, "f": f, "g": g}

          # Create an NLP solver
          # Function solver = nlpsol("solver", "sqpmethod", nlp, opts)
          solver = nlpsol("solver", "feasiblesqpmethod", nlp, opts)

          # Solve the problem
          arg = {}
          arg["x0"] = 1.0 #[2.5,3.0,0.75}
          arg["lbg"] = -5
          arg["ubg"] = 5
          res = solver(**arg)
          self.assertEqual(solver.stats()["iter_count"],iter_count)
          self.checkarray(res["x"],[-4],digits=6)

if __name__ == '__main__':
    unittest.main()
