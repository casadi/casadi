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


      #include "bspline_interpolant.hpp"
      #include <string>

      const std::string casadi::BSplineInterpolant::meta_doc=
      "\n"
"\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_239 \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------------+--------------+-----------------------------------+\n"
"|          Id           |     Type     |            Description            |\n"
"+=======================+==============+===================================+\n"
"| algorithm             | OT_STRING    | Algorithm used for fitting the    |\n"
"|                       |              | data: 'not_a_knot' (default, same |\n"
"|                       |              | as Matlab), 'smooth_linear'.      |\n"
"+-----------------------+--------------+-----------------------------------+\n"
"| degree                | OT_INTVECTOR | Sets, for each grid dimension,    |\n"
"|                       |              | the degree of the spline.         |\n"
"+-----------------------+--------------+-----------------------------------+\n"
"| linear_solver         | OT_STRING    | Solver used for constructing the  |\n"
"|                       |              | coefficient tensor.               |\n"
"+-----------------------+--------------+-----------------------------------+\n"
"| linear_solver_options | OT_DICT      | Options to be passed to the       |\n"
"|                       |              | linear solver.                    |\n"
"+-----------------------+--------------+-----------------------------------+\n"
"| smooth_linear_frac    | OT_DOUBLE    | When 'smooth_linear' algorithm is |\n"
"|                       |              | active, determines sharpness      |\n"
"|                       |              | between 0 (sharp, as linear       |\n"
"|                       |              | interpolation) and 0.5            |\n"
"|                       |              | (smooth).Default value is 0.1.    |\n"
"+-----------------------+--------------+-----------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
