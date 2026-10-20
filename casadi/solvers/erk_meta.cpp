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


      #include "erk.hpp"
      #include <string>

      const std::string casadi::Erk::meta_doc=
      "\n"
"\n"
"\n"
"Explicit Runge-Kutta integrator with fixed steps for ODEs: the classical \n"
"method of order 4.\n"
"\n"
"Steps are independent of the output times: each run of steps over a \n"
"control group (or up to the end) is divided into equal steps \n"
"of at most step_size, and outputs inside a step come from a dense output \n"
"of order 3 at no extra cost. Forward sensitivities are the exact \n"
"derivatives of the discrete solution and adjoint sensitivities are the \n"
"exact discrete adjoint of the taped steps. The tape has a known bound, so \n"
"adjoints can be code-generated.\n"
"\n"
"\n"
;
