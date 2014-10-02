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


      #include "snopt_interface.hpp"
      #include <string>

      const std::string casadi::SnoptInterface::meta_doc=
      "\n"
"SNOPT interface\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"|       Id        |      Type       |     Default     |   Description   |\n"
"+=================+=================+=================+=================+\n"
"| detect_linear   | OT_BOOLEAN      | true            | Make an effort  |\n"
"|                 |                 |                 | to treat linear |\n"
"|                 |                 |                 | constraints and |\n"
"|                 |                 |                 | linear          |\n"
"|                 |                 |                 | variables       |\n"
"|                 |                 |                 | specially.      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| print file      | OT_STRING       |                 |                 |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| print_time      | OT_BOOLEAN      | true            | print           |\n"
"|                 |                 |                 | information     |\n"
"|                 |                 |                 | about execution |\n"
"|                 |                 |                 | time            |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| specs file      | OT_STRING       |                 |                 |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| start           | OT_STRING       | \"Cold\"          | (Cold|Basis|War |\n"
"|                 |                 |                 | m)              |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| summary         | OT_BOOLEAN      | true            |                 |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"\n"
"\n"
">List of available monitors\n"
"\n"
"+-----------+\n"
"|    Id     |\n"
"+===========+\n"
"| eval_nlp  |\n"
"+-----------+\n"
"| setup_nlp |\n"
"+-----------+\n"
"\n"
"\n"
">List of available stats\n"
"\n"
"+----------------+\n"
"|       Id       |\n"
"+================+\n"
"| iter_count     |\n"
"+----------------+\n"
"| iterations     |\n"
"+----------------+\n"
"| n_callback_fun |\n"
"+----------------+\n"
"| n_eval_grad_f  |\n"
"+----------------+\n"
"| n_eval_jac_g   |\n"
"+----------------+\n"
"| return_status  |\n"
"+----------------+\n"
"| t_callback_fun |\n"
"+----------------+\n"
"| t_eval_grad_f  |\n"
"+----------------+\n"
"| t_eval_jac_g   |\n"
"+----------------+\n"
"| t_mainloop     |\n"
"+----------------+\n"
"\n"
"\n"
"\n"
"\n"
;
