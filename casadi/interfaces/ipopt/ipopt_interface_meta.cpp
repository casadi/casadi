/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


      #include "ipopt_interface.hpp"
      #include <string>

      const std::string casadi::IpoptInterface::meta_doc=
      "\n"
"\n"
"\n"
"When in warmstart mode, output NLPSOL_LAM_X may be used as input\n"
"\n"
"NOTE: Even when max_iter == 0, it is not guaranteed that \n"
"input(NLPSOL_X0) == output(NLPSOL_X). Indeed if bounds on X or \n"
"constraints are unmet, they will differ.\n"
"\n"
"For a good tutorial on IPOPT, see http://drops.dagstuhl.de/volltexte/2009/2089/pdf/09061.WaechterAndreas.Paper.2089.pdf \n"
"\n"
"A good resource about the algorithms in IPOPT is: Wachter and L. T. \n"
"Biegler, On the Implementation of an Interior-Point Filter Line-Search\n"
" Algorithm for Large-Scale Nonlinear Programming, Mathematical \n"
"Programming 106(1), pp. 25-57, 2006 (As Research Report RC 23149, IBM \n"
"T. J. Watson Research Center, Yorktown, USA\n"
"\n"
"Caveats:\n"
"with default options, multipliers for the decision variables are wrong\n"
" for equality constraints. Change the 'fixed_variable_treatment' to \n"
"'make_constraint' or 'relax_bounds' to obtain correct results.\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_21y \n"
"\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+--------------------------+-------------+---------------------------------+\n"
"|            Id            |    Type     |           Description           |\n"
"+==========================+=============+=================================+\n"
"| clip_inactive_lam        | OT_BOOL     | Explicitly set Lagrange         |\n"
"|                          |             | multipliers to 0 when bound is  |\n"
"|                          |             | deemed inactive (default:       |\n"
"|                          |             | false).                         |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| con_integer_md           | OT_DICT     | Integer metadata (a dictionary  |\n"
"|                          |             | with lists of integers) about   |\n"
"|                          |             | constraints to be passed to     |\n"
"|                          |             | IPOPT                           |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| con_numeric_md           | OT_DICT     | Numeric metadata (a dictionary  |\n"
"|                          |             | with lists of reals) about      |\n"
"|                          |             | constraints to be passed to     |\n"
"|                          |             | IPOPT                           |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| con_string_md            | OT_DICT     | String metadata (a dictionary   |\n"
"|                          |             | with lists of strings) about    |\n"
"|                          |             | constraints to be passed to     |\n"
"|                          |             | IPOPT                           |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| convexify_margin         | OT_DOUBLE   | When using a convexification    |\n"
"|                          |             | strategy, make sure that the    |\n"
"|                          |             | smallest eigenvalue is at least |\n"
"|                          |             | this (default: 1e-7).           |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| convexify_strategy       | OT_STRING   | NONE|regularize|eigen-          |\n"
"|                          |             | reflect|eigen-clip. Strategy to |\n"
"|                          |             | convexify the Lagrange Hessian  |\n"
"|                          |             | before passing it to the        |\n"
"|                          |             | solver.                         |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| grad_f                   | OT_FUNCTION | Function for calculating the    |\n"
"|                          |             | gradient of the objective       |\n"
"|                          |             | (column, autogenerated by       |\n"
"|                          |             | default)                        |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| hess_lag                 | OT_FUNCTION | Function for calculating the    |\n"
"|                          |             | Hessian of the Lagrangian       |\n"
"|                          |             | (autogenerated by default)      |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| inactive_lam_strategy    | OT_STRING   | Strategy to detect if a bound   |\n"
"|                          |             | is inactive. RELTOL: use        |\n"
"|                          |             | solver-defined constraint       |\n"
"|                          |             | tolerance *                     |\n"
"|                          |             | inactive_lam_value|abstol: use  |\n"
"|                          |             | inactive_lam_value              |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| inactive_lam_value       | OT_DOUBLE   | Value used in                   |\n"
"|                          |             | inactive_lam_strategy (default: |\n"
"|                          |             | 10).                            |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| ipopt                    | OT_DICT     | Options to be passed to IPOPT   |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| jac_g                    | OT_FUNCTION | Function for calculating the    |\n"
"|                          |             | Jacobian of the constraints     |\n"
"|                          |             | (autogenerated by default)      |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| max_iter_eig             | OT_DOUBLE   | Maximum number of iterations to |\n"
"|                          |             | compute an eigenvalue           |\n"
"|                          |             | decomposition (default: 50).    |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| pass_nonlinear_variables | OT_BOOL     | Pass list of variables entering |\n"
"|                          |             | nonlinearly to IPOPT            |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| var_integer_md           | OT_DICT     | Integer metadata (a dictionary  |\n"
"|                          |             | with lists of integers) about   |\n"
"|                          |             | variables to be passed to IPOPT |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| var_numeric_md           | OT_DICT     | Numeric metadata (a dictionary  |\n"
"|                          |             | with lists of reals) about      |\n"
"|                          |             | variables to be passed to IPOPT |\n"
"+--------------------------+-------------+---------------------------------+\n"
"| var_string_md            | OT_DICT     | String metadata (a dictionary   |\n"
"|                          |             | with lists of strings) about    |\n"
"|                          |             | variables to be passed to IPOPT |\n"
"+--------------------------+-------------+---------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
