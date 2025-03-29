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


      #include "slicot_dple.hpp"
      #include <string>

      const std::string casadi::SlicotDple::meta_doc=
      "\n"
"\n"
"\n"
"An efficient solver for Discrete Periodic Lyapunov Equations using \n"
"SLICOT\n"
"\n"
"Uses Periodic Schur Decomposition ('psd') and does not assume positive\n"
" definiteness. Based on Periodic Lyapunov equations: some applications\n"
" and new algorithms. Int. J. Control, vol. 67, pp. 69-87, 1997.\n"
"\n"
"Overview of the method: J. Gillis Practical Methods for Approximate \n"
"Robust Periodic Optimal Control ofNonlinear Mechanical Systems, PhD \n"
"Thesis, KULeuven, 2015\n"
"\n"
"Extra doc: https://github.com/casadi/casadi/wiki/L_22j \n"
"\n"
"\n"
">List of available options\n"
"\n"
"+-----------------------+-----------+--------------------------------------+\n"
"|          Id           |   Type    |             Description              |\n"
"+=======================+===========+======================================+\n"
"| linear_solver         | OT_STRING | User-defined linear solver class.    |\n"
"|                       |           | Needed for sensitivities.            |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| linear_solver_options | OT_DICT   | Options to be passed to the linear   |\n"
"|                       |           | solver.                              |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"| psd_num_zero          | OT_DOUBLE | Numerical zero used in Periodic      |\n"
"|                       |           | Schur decomposition with slicot.This |\n"
"|                       |           | option is needed when your systems   |\n"
"|                       |           | has Floquet multiplierszero or close |\n"
"|                       |           | to zero                              |\n"
"+-----------------------+-----------+--------------------------------------+\n"
"\n"
"\n"
"\n"
"\n"
;
