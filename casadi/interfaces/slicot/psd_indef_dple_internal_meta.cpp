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


      #include "psd_indef_dple_internal.hpp"
      #include <string>

      const std::string casadi::PsdIndefDpleInternal::meta_doc=
      "\n"
"An efficient solver for Discrete Periodic Lyapunov Equations using\n"
"SLICOT\n"
"\n"
"Uses Periodic Schur Decomposition ('psd') and does not assume positive\n"
"definiteness. Based on Periodic Lyapunov equations: some applications\n"
"and new algorithms. Int. J. Control, vol. 67, pp. 69-87, 1997.\n"
"\n"
"\n"
;
