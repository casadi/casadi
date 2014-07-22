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


      #include "ooqp_internal.hpp"
      #include <string>

      const std::string casadi::OOQPInternal::meta_doc=
      "\n"
"\n"
">List of available options\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"|       Id        |      Type       |     Default     |   Description   |\n"
"+=================+=================+=================+=================+\n"
"| ad_mode         | OT_STRING       | \"automatic\"     | How to          |\n"
"|                 |                 |                 | calculate the   |\n"
"|                 |                 |                 | Jacobians.      |\n"
"|                 |                 |                 | (forward: only  |\n"
"|                 |                 |                 | forward         |\n"
"|                 |                 |                 | mode|reverse:   |\n"
"|                 |                 |                 | only adjoint    |\n"
"|                 |                 |                 | mode|automatic: |\n"
"|                 |                 |                 | a heuristic     |\n"
"|                 |                 |                 | decides which   |\n"
"|                 |                 |                 | is more         |\n"
"|                 |                 |                 | appropriate)    |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| artol           | OT_REAL         | 0.000           | tolerance as    |\n"
"|                 |                 |                 | provided with   |\n"
"|                 |                 |                 | setArTol to     |\n"
"|                 |                 |                 | OOQP            |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| derivative_gene | OT_DERIVATIVEGE | GenericType()   | Function that   |\n"
"| rator           | NERATOR         |                 | returns a       |\n"
"|                 |                 |                 | derivative      |\n"
"|                 |                 |                 | function given  |\n"
"|                 |                 |                 | a number of     |\n"
"|                 |                 |                 | forward and     |\n"
"|                 |                 |                 | reverse         |\n"
"|                 |                 |                 | directional     |\n"
"|                 |                 |                 | derivative,     |\n"
"|                 |                 |                 | overrides       |\n"
"|                 |                 |                 | internal        |\n"
"|                 |                 |                 | routines. Check |\n"
"|                 |                 |                 | documentation   |\n"
"|                 |                 |                 | of DerivativeGe |\n"
"|                 |                 |                 | nerator.        |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| gather_stats    | OT_BOOLEAN      | false           | Flag to         |\n"
"|                 |                 |                 | indicate        |\n"
"|                 |                 |                 | whether         |\n"
"|                 |                 |                 | statistics must |\n"
"|                 |                 |                 | be gathered     |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| inputs_check    | OT_BOOLEAN      | true            | Throw           |\n"
"|                 |                 |                 | exceptions when |\n"
"|                 |                 |                 | the numerical   |\n"
"|                 |                 |                 | values of the   |\n"
"|                 |                 |                 | inputs don't    |\n"
"|                 |                 |                 | make sense      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| monitor         | OT_STRINGVECTOR | GenericType()   | Monitors to be  |\n"
"|                 |                 |                 | activated (inpu |\n"
"|                 |                 |                 | ts|outputs)     |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| mutol           | OT_REAL         | 0.000           | tolerance as    |\n"
"|                 |                 |                 | provided with   |\n"
"|                 |                 |                 | setMuTol to     |\n"
"|                 |                 |                 | OOQP            |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| name            | OT_STRING       | \"unnamed_shared | name of the     |\n"
"|                 |                 | _object\"        | object          |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| print_level     | OT_INTEGER      | 0               | Print level.    |\n"
"|                 |                 |                 | OOQP listens to |\n"
"|                 |                 |                 | print_level 0,  |\n"
"|                 |                 |                 | 10 and 100      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| regularity_chec | OT_BOOLEAN      | true            | Throw           |\n"
"| k               |                 |                 | exceptions when |\n"
"|                 |                 |                 | NaN or Inf      |\n"
"|                 |                 |                 | appears during  |\n"
"|                 |                 |                 | evaluation      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| user_data       | OT_VOIDPTR      | GenericType()   | A user-defined  |\n"
"|                 |                 |                 | field that can  |\n"
"|                 |                 |                 | be used to      |\n"
"|                 |                 |                 | identify the    |\n"
"|                 |                 |                 | function or     |\n"
"|                 |                 |                 | pass additional |\n"
"|                 |                 |                 | information     |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| verbose         | OT_BOOLEAN      | false           | Verbose         |\n"
"|                 |                 |                 | evaluation  for |\n"
"|                 |                 |                 | debugging       |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"\n"
"\n"
"\n"
"\n"
;
