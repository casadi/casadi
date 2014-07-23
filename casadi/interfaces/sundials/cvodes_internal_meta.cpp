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


      #include "cvodes_internal.hpp"
      #include <string>

      const std::string casadi::CVodesInternal::meta_doc=
      "\n"
"Interface to CVodes from the Sundials suite.\n"
"\n"
"A call to evaluate will integrate to the end.\n"
"\n"
"You can retrieve the entire state trajectory as follows, after the\n"
"evaluate call: Call reset. Then call integrate(t_i) and getOuput for a\n"
"series of times t_i.\n"
"\n"
"\n"
;
