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


      #include "ipmc_interface.hpp"
      #include <string>

      const std::string casadi::IpmcInterface::meta_doc=
      "\n"
"\n"
">List of available options\n"
"\n"
"+---------------------+--------------+-------------------------------------+\n"
"|         Id          |     Type     |             Description             |\n"
"+=====================+==============+=====================================+\n"
"| N                   | OT_INT       | OCP horizon                         |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| debug               | OT_BOOL      | Write the expected and actual       |\n"
"|                     |              | Jacobian structure to               |\n"
"|                     |              | debug_ipmc_*.mtx [false].           |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| ipmc                | OT_DICT      | Options to be passed to ipmc.       |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| ng                  | OT_INTVECTOR | Number of non-dynamic constraints,  |\n"
"|                     |              | length N+1                          |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| nu                  | OT_INTVECTOR | Number of controls, length N+1      |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| nx                  | OT_INTVECTOR | Number of states, length N+1        |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| nxc                 | OT_INT       | Number of trailing states per stage |\n"
"|                     |              | that are constant, x_{k+1}=x_k;     |\n"
"|                     |              | ipmc exploits them in the Riccati   |\n"
"|                     |              | recursion. Needs                    |\n"
"|                     |              | structure_detection 'manual' or     |\n"
"|                     |              | 'auto' [0].                         |\n"
"+---------------------+--------------+-------------------------------------+\n"
"| structure_detection | OT_STRING    | Structure detection: none, auto or  |\n"
"|                     |              | manual [none].                      |\n"
"+---------------------+--------------+-------------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
