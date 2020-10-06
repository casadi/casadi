/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
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

namespace {
const std::string meta_doc_0 =
"\n"
"When in warmstart mode, output NLPSOL_LAM_X may be used as input\n"
"\n"
"NOTE: Even when max_iter == 0, it is not guaranteed that\n"
"input(NLPSOL_X0) == output(NLPSOL_X). Indeed if bounds on X or\n"
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
">List of available options\n"
"\n"
"Too long list to be displayed."
;
const std::string meta_doc_1 =
"Too long list to be displayed."
;
const std::string meta_doc_2 =
"Too long list to be displayed."
"\n"
">List of available monitors\n"
"\n"
"+-------------+\n"
"|     Id      |\n"
"+=============+\n"
"| eval_f      |\n"
"+-------------+\n"
"| eval_g      |\n"
"+-------------+\n"
"| eval_grad_f |\n"
"+-------------+\n"
"| eval_h      |\n"
"+-------------+\n"
"| eval_jac_g  |\n"
"+-------------+\n"
"\n"
"\n"
">List of available stats\n"
"\n"
"+--------------------+\n"
"|         Id         |\n"
"+====================+\n"
"| con_integer_md     |\n"
"+--------------------+\n"
"| con_numeric_md     |\n"
"+--------------------+\n"
"| con_string_md      |\n"
"+--------------------+\n"
"| iter_count         |\n"
"+--------------------+\n"
"| iteration          |\n"
"+--------------------+\n"
"| iterations         |\n"
"+--------------------+\n"
"| n_eval_f           |\n"
"+--------------------+\n"
"| n_eval_g           |\n"
"+--------------------+\n"
"| n_eval_grad_f      |\n"
"+--------------------+\n"
"| n_eval_h           |\n"
"+--------------------+\n"
"| n_eval_jac_g       |\n"
"+--------------------+\n"
"| return_status      |\n"
"+--------------------+\n"
"| t_callback_fun     |\n"
"+--------------------+\n"
"| t_callback_prepare |\n"
"+--------------------+\n"
"| t_eval_f           |\n"
"+--------------------+\n"
"| t_eval_g           |\n"
"+--------------------+\n"
"| t_eval_grad_f      |\n"
"+--------------------+\n"
"| t_eval_h           |\n"
"+--------------------+\n"
"| t_eval_jac_g       |\n"
"+--------------------+\n"
"| t_mainloop         |\n"
"+--------------------+\n"
"| var_integer_md     |\n"
"+--------------------+\n"
"| var_numeric_md     |\n"
"+--------------------+\n"
"| var_string_md      |\n"
"+--------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
}

const std::string casadi::IpoptInterface::meta_doc = meta_doc_0 + meta_doc_1 + meta_doc_2;
