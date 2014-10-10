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
# -*- coding: utf-8 -*-
from casadi import *
import numpy as NP
import matplotlib.pyplot as plt

""" 
Example program demonstrating sensitivity analysis with sIPOPT from CasADi

Test problem (Ganesh & Biegler, A reduced Hessian strategy for sensitivity analysis of optimal flowsheets, AIChE 33, 1987, pp. 282-296)

  min     x1^2 + x2^2 + x3^2
  s.t.    6*x1 + 3&x2 + 2*x3 - pi = 0
          p2*x1 + x2 - x3 - 1 = 0
          x1, x2, x3 >= 0

@author Joel Andersson, K.U. Leuven 2012
"""

# Optimization variables
x = SX.sym("x",3)
  
 # Parameters
p = SX.sym("p",2)
  
# Objective
f = x[0]*x[0] + x[1]*x[1] + x[2]*x[2]
  
# Constraints
g = vertcat(( \
       6*x[0] + 3*x[1] + 2*x[2] - p[0],
    p[1]*x[0] +   x[1] -   x[2] -    1))
  
# Augment the parameters to the list of variables
x = vertcat((x,p))
  
# Fix the parameters by additional equations
g = vertcat((g,p))
  
# Original parameter values
p_a  = [5.00,1.00]
  
# Perturbed parameter values
p_b  = [4.50,1.00]

# Initial guess and bounds for the optimization variables
x0  = [0.15, 0.15, 0.00, p_a[0], p_a[1]]
lbx = [0.00, 0.00, 0.00,   -inf,   -inf]
ubx = [ inf,  inf,  inf,    inf,    inf]
  
# Nonlinear bounds
lbg = [0.00, 0.00, p_a[0], p_a[1]]
ubg = [0.00, 0.00, p_a[0], p_a[1]]
    
# Create NLP solver
nlp = SXFunction(nlpIn(x=x),nlpOut(f=f,g=g))
solver = NlpSolver("ipopt", nlp)
  
# Mark the parameters amongst the constraints (see sIPOPT documentation)
con_integer_md = {}
con_integer_md["sens_init_constr"] = [0,0,1,2]
solver.setOption("con_integer_md",con_integer_md)
  
# Mark the parameters amongst the variables (see sIPOPT documentation)
var_integer_md = {}
var_integer_md["sens_state_1"] = [0,0,0,1,2]
solver.setOption("var_integer_md",var_integer_md)

# Pass the perturbed values (see sIPOPT documentation)
var_numeric_md = {}
var_numeric_md["sens_state_value_1"] = [0,0,0,p_b[0],p_b[1]]
solver.setOption("var_numeric_md",var_numeric_md)
  
# Enable sensitivities
solver.setOption("run_sens","yes")
solver.setOption("n_sens_steps", 1)
  
# Initialize solver
solver.init()
  
# Solve NLP
solver.setInput( x0, "x0")
solver.setInput(lbx, "lbx")
solver.setInput(ubx, "ubx")
solver.setInput(lbg, "lbg")
solver.setInput(ubg, "ubg")
solver.evaluate()
  
# Print the solution
print "----"
print "Minimal cost " , solver.getOutput("f")
print "----"

print "Nominal solution"
print "x = " , solver.output("x").data()
print "----"
  
print "perturbed solution"
var_numeric_md = solver.getStat("var_numeric_md")
print "x = " , var_numeric_md["sens_sol_state_1"]
print "----"
  
print "Dual bound multipliers"
print "z_L = " , var_numeric_md["sens_sol_state_1_z_L"]
print "z_U = " , var_numeric_md["sens_sol_state_1_z_U"]
print "----"
  
print "Constraint multipliers"
con_numeric_md = solver.getStat("con_numeric_md")
print "lambda = " , con_numeric_md["sens_sol_state_1"]
print "----"
  
