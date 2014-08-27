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


      #include "condensing_indef_dple_internal.hpp"
      #include <string>

      const std::string casadi::CondensingIndefDpleInternal::meta_doc=
      "\n"
"Solving the Discrete Periodic Lyapunov Equations by condensing the\n"
"entire period to a single Discrete Lyapunov Equation\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"|       Id        |      Type       |     Default     |   Description   |\n"
"+=================+=================+=================+=================+\n"
"| dle_solver      | OT_STRING       | GenericType()   | User-defined    |\n"
"|                 |                 |                 | Dle solver      |\n"
"|                 |                 |                 | class. Needed   |\n"
"|                 |                 |                 | for             |\n"
"|                 |                 |                 | sensitivities.  |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| dle_solver_opti | OT_DICTIONARY   | GenericType()   | Options to be   |\n"
"| ons             |                 |                 | passed to the   |\n"
"|                 |                 |                 | Dle solver.     |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"\n"
"\n"
"\n"
"\n"
;
