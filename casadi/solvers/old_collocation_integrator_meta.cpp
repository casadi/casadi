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


      #include "old_collocation_integrator.hpp"
      #include <string>

      const std::string casadi::OldCollocationIntegrator::meta_doc=
      "\n"
"Collocation integrator ODE/DAE integrator based on collocation\n"
"\n"
"The method is still under development\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"|       Id        |      Type       |     Default     |   Description   |\n"
"+=================+=================+=================+=================+\n"
"| collocation_sch | OT_STRING       | \"radau\"         | Collocation     |\n"
"| eme             |                 |                 | scheme (radau|l |\n"
"|                 |                 |                 | egendre)        |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| expand_f        | OT_BOOLEAN      | false           | Expand the      |\n"
"|                 |                 |                 | ODE/DAE         |\n"
"|                 |                 |                 | residual        |\n"
"|                 |                 |                 | function in an  |\n"
"|                 |                 |                 | SX graph        |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| expand_q        | OT_BOOLEAN      | false           | Expand the      |\n"
"|                 |                 |                 | quadrature      |\n"
"|                 |                 |                 | function in an  |\n"
"|                 |                 |                 | SX graph        |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| hotstart        | OT_BOOLEAN      | true            | Initialize the  |\n"
"|                 |                 |                 | trajectory at   |\n"
"|                 |                 |                 | the previous    |\n"
"|                 |                 |                 | solution        |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| implicit_solver | OT_STRING       | GenericType()   | An implicit     |\n"
"|                 |                 |                 | function solver |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| implicit_solver | OT_DICTIONARY   | GenericType()   | Options to be   |\n"
"| _options        |                 |                 | passed to the   |\n"
"|                 |                 |                 | implicit solver |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| interpolation_o | OT_INTEGER      | 3               | Order of the    |\n"
"| rder            |                 |                 | interpolating   |\n"
"|                 |                 |                 | polynomials     |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| number_of_finit | OT_INTEGER      | 20              | Number of       |\n"
"| e_elements      |                 |                 | finite elements |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| startup_integra | OT_STRING       | GenericType()   | An ODE/DAE      |\n"
"| tor             |                 |                 | integrator that |\n"
"|                 |                 |                 | can be used to  |\n"
"|                 |                 |                 | generate a      |\n"
"|                 |                 |                 | startup         |\n"
"|                 |                 |                 | trajectory      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"| startup_integra | OT_DICTIONARY   | GenericType()   | Options to be   |\n"
"| tor_options     |                 |                 | passed to the   |\n"
"|                 |                 |                 | startup         |\n"
"|                 |                 |                 | integrator      |\n"
"+-----------------+-----------------+-----------------+-----------------+\n"
"\n"
"\n"
"\n"
"\n"
;
