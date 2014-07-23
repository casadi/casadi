/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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


      #include "ipopt_internal.hpp"
      #include <string>

      const std::string casadi::IpoptInternal::meta_doc=
      "\n"
"When in warmstart mode, output NLP_SOLVER_LAM_X may be used as input\n"
"\n"
"NOTE: Even when max_iter == 0, it is not guaranteed that\n"
"input(NLP_SOLVER_X0) == output(NLP_SOLVER_X). Indeed if bounds on X or\n"
"constraints are unmet, they will differ.\n"
"\n"
"For a good tutorial on IPOPT, seehttp://drops.dagstuhl.de/volltexte/2009/2089/pdf/09061.WaechterAndreas.Paper.2089.pdf\n"
"\n"
"A good resource about the algorithms in IPOPT is: Wachter and L. T.\n"
"Biegler, On the Implementation of an Interior-Point Filter Line-Search\n"
"Algorithm for Large-Scale Nonlinear Programming, Mathematical\n"
"Programming 106(1), pp. 25-57, 2006 (As Research Report RC 23149, IBM\n"
"T. J. Watson Research Center, Yorktown, USA\n"
"\n"
"Caveats:\n"
"with default options, multipliers for the decision variables are wrong\n"
"for equality constraints. Change the 'fixed_variable_treatment' to\n"
"'make_constraint' or 'relax_bounds' to obtain correct results.\n"
"\n"
"\n"
"\n"
;
