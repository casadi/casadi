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


      #include "collocation.hpp"
      #include <string>

      const std::string casadi::Collocation::meta_doc=
      "\n"
"\n"
"\n"
"Fixed-step implicit Runge-Kutta integrator ODE/DAE integrator based \n"
"on collocation schemes\n"
"\n"
"The method is still under development\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_234 \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+---------------------------+-----------+----------------------------------+\n"
"|            Id             |   Type    |           Description            |\n"
"+===========================+===========+==================================+\n"
"| collocation_scheme        | OT_STRING | Collocation scheme:              |\n"
"|                           |           | radau|legendre                   |\n"
"+---------------------------+-----------+----------------------------------+\n"
"| interpolation_order       | OT_INT    | Order of the interpolating       |\n"
"|                           |           | polynomials                      |\n"
"+---------------------------+-----------+----------------------------------+\n"
"| number_of_finite_elements | OT_INT    | Target number of finite          |\n"
"|                           |           | elements. The actual number may  |\n"
"|                           |           | be higher to accommodate all     |\n"
"|                           |           | output times                     |\n"
"+---------------------------+-----------+----------------------------------+\n"
"| rootfinder                | OT_STRING | An implicit function solver      |\n"
"+---------------------------+-----------+----------------------------------+\n"
"| rootfinder_options        | OT_DICT   | Options to be passed to the NLP  |\n"
"|                           |           | Solver                           |\n"
"+---------------------------+-----------+----------------------------------+\n"
"| simplify                  | OT_BOOL   | Implement as MX Function         |\n"
"|                           |           | (codegeneratable/serializable)   |\n"
"|                           |           | default: false                   |\n"
"+---------------------------+-----------+----------------------------------+\n"
"| simplify_options          | OT_DICT   | Any options to pass to           |\n"
"|                           |           | simplified form Function         |\n"
"|                           |           | constructor                      |\n"
"+---------------------------+-----------+----------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
