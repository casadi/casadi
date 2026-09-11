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


      #include "shf_interpolant.hpp"
      #include <string>

      const std::string casadi::SHFInterpolant::meta_doc=
      "\n"
"\n"
"\n"
"Smooth hat function (SHF) approximation of the piecewise-linear interpolant: C^k, \n"
"coincides with the look-up table outside epsilon-balls around the interior grid \n"
"points, no fitting step. Wraps shf_spline.\n"
"\n"
"\n"
">List of available options\n"
"\n"
"+--------------------+----------------+-------------------------------------------+\n"
"|         Id         |      Type      |                Description                |\n"
"+====================+================+===========================================+\n"
"| epsilon            | OT_DOUBLEVECTOR| Absolute half-width of the smoothing      |\n"
"|                    |                | balls around the interior grid points,    |\n"
"|                    |                | one value shared by all axes or one per   |\n"
"|                    |                | axis. Default: 0.1 times the smallest     |\n"
"|                    |                | grid spacing of each axis.                |\n"
"+--------------------+----------------+-------------------------------------------+\n"
"| epsilon_parametric | OT_BOOL        | Take epsilon as an additional (non-       |\n"
"|                    |                | differentiable) input 'eps' of length     |\n"
"|                    |                | ndim instead of fixing it at              |\n"
"|                    |                | construction. Default: false.             |\n"
"+--------------------+----------------+-------------------------------------------+\n"
"| smoothness_order   | OT_INT         | Smoothness order k: the approximation is  |\n"
"|                    |                | C^k. Default: 3.                          |\n"
"+--------------------+----------------+-------------------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
