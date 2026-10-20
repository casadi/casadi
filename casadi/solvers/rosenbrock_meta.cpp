/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl,
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


      #include "rosenbrock.hpp"
      #include <string>

      const std::string casadi::Rosenbrock::meta_doc=
      "\n"
"\n"
"\n"
"Adaptive Rosenbrock-Wanner integrator for stiff ODEs and semi-explicit \n"
"index-1 DAEs.\n"
"\n"
"An 8-stage, stiffly accurate, L-stable method of order 5 with an \n"
"embedded method of order 4. Each step needs one Jacobian and one sparse \n"
"LU factorization; there is no Newton iteration. Step sizes are chosen \n"
"from the nominal solution only, so forward sensitivities are the exact \n"
"derivatives of the discrete solution (internal numerical \n"
"differentiation) and adjoint sensitivities are the exact discrete \n"
"adjoint of the taped steps.\n"
"\n"
"\n"
;
