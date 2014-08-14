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


      #include "qpoases_interface.hpp"
      #include <string>

      const std::string casadi::QpoasesInterface::meta_doc=
      "\n"
"Interface to QPOases Solver for quadratic programming\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"|       Id        |      Type       |     Default     |   Description   |\n"
"+=================+=================+=================+=================+\n"
"| CPUtime         | OT_REAL         | GenericType()   | The maximum     |\n"
"|                 |                 |                 | allowed CPU     |\n"
"|                 |                 |                 | time in seconds |\n"
"|                 |                 |                 | for the whole   |\n"
"|                 |                 |                 | initialisation  |\n"
"|                 |                 |                 | (and the        |\n"
"|                 |                 |                 | actually        |\n"
"|                 |                 |                 | required one on |\n"
"|                 |                 |                 | output).        |\n"
"|                 |                 |                 | Disabled if     |\n"
"|                 |                 |                 | unset.          |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| enableEqualitie | OT_BOOLEAN      | BooleanType_to_ | Specifies       |\n"
"| s               |                 | bool            | whether         |\n"
"|                 |                 |                 | equalities      |\n"
"|                 |                 |                 | should be       |\n"
"|                 |                 |                 | treated as      |\n"
"|                 |                 |                 | always active   |\n"
"|                 |                 |                 | (True) or not   |\n"
"|                 |                 |                 | (False)         |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| enableFarBounds | OT_BOOLEAN      | BooleanType_to_ | Enables the use |\n"
"|                 |                 | bool            | of far bounds.  |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| enableFlippingB | OT_BOOLEAN      | BooleanType_to_ | Enables the use |\n"
"| ounds           |                 | bool            | of flipping     |\n"
"|                 |                 |                 | bounds.         |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| enableFullLITes | OT_BOOLEAN      | BooleanType_to_ | Enables         |\n"
"| ts              |                 | bool            | condition-      |\n"
"|                 |                 |                 | hardened (but   |\n"
"|                 |                 |                 | more expensive) |\n"
"|                 |                 |                 | LI test.        |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| enableNZCTests  | OT_BOOLEAN      | BooleanType_to_ | Enables nonzero |\n"
"|                 |                 | bool            | curvature       |\n"
"|                 |                 |                 | tests.          |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| enableRamping   | OT_BOOLEAN      | BooleanType_to_ | Enables         |\n"
"|                 |                 | bool            | ramping.        |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| enableRegularis | OT_BOOLEAN      | BooleanType_to_ | Enables         |\n"
"| ation           |                 | bool            | automatic       |\n"
"|                 |                 |                 | Hessian         |\n"
"|                 |                 |                 | regularisation. |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| initialStatusBo | OT_STRING       | SubjectToStatus | Initial status  |\n"
"| unds            |                 | _to_string      | of bounds at    |\n"
"|                 |                 |                 | first           |\n"
"|                 |                 |                 | iteration.      |\n"
"|                 |                 |                 | (inactive::all  |\n"
"|                 |                 |                 | bounds inactive |\n"
"|                 |                 |                 | |lower::all     |\n"
"|                 |                 |                 | bounds active   |\n"
"|                 |                 |                 | at their lower  |\n"
"|                 |                 |                 | bound|upper::al |\n"
"|                 |                 |                 | l bounds active |\n"
"|                 |                 |                 | at their upper  |\n"
"|                 |                 |                 | bound)          |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| nWSR            | OT_INTEGER      | GenericType()   | The maximum     |\n"
"|                 |                 |                 | number of       |\n"
"|                 |                 |                 | working set     |\n"
"|                 |                 |                 | recalculations  |\n"
"|                 |                 |                 | to be performed |\n"
"|                 |                 |                 | during the      |\n"
"|                 |                 |                 | initial         |\n"
"|                 |                 |                 | homotopy.       |\n"
"|                 |                 |                 | Default is 5(nx |\n"
"|                 |                 |                 | + nc)           |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| printLevel      | OT_STRING       | PrintLevel_to_s | Defines the     |\n"
"|                 |                 | tring           | amount of text  |\n"
"|                 |                 |                 | output during   |\n"
"|                 |                 |                 | QP solution,    |\n"
"|                 |                 |                 | see Section 5.7 |\n"
"|                 |                 |                 | (none|low|mediu |\n"
"|                 |                 |                 | m|high)         |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"\n"
"\n"
"\n"
"\n"
;
