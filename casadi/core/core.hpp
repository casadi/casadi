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


#ifndef CASADI_CORE_HPP
#define CASADI_CORE_HPP

// Scalar expressions (why do I need to put it up here?)
#include "sx/sx_elem.hpp"

// Generic tools
#include "polynomial.hpp"
#include "std_vector_tools.hpp"
#include "global_options.hpp"
#include "casadi_meta.hpp"

// Matrices
#include "matrix.hpp"

// Matrix expressions
#include "mx/mx.hpp"

// Functions
#include "function/oracle.hpp"
#include "function/code_generator.hpp"
#include "function/compiler.hpp"
#include "function/callback.hpp"
#include "function/integrator.hpp"
#include "function/qpsol.hpp"
#include "function/nlpsol.hpp"
#include "function/rootfinder.hpp"
#include "function/linsol.hpp"
#include "function/jit.hpp"
#include "function/external.hpp"

// Misc
#include "misc/integration_tools.hpp"
#include "misc/nlp_builder.hpp"
#include "misc/variable.hpp"
#include "misc/dae_builder.hpp"
#include "misc/xml_file.hpp"

#endif // CASADI_CORE_HPP
