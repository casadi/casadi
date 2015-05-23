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


#include "sx_tools.hpp"
#include "../function/sx_function_internal.hpp"
#include "../casadi_math.hpp"
#include "../casadi_exception.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../std_vector_tools.hpp"
using namespace std;

namespace casadi {

  Matrix<double> evalf(const SX &ex, const SX &v, const Matrix<double> &vdef) {
    SXFunction fcn("tmp", make_vector(v), make_vector(ex));
    fcn.setInput(vdef);
    fcn.evaluate();
    return fcn.output();
  }

  Matrix<double> evalf(const SX& ex) {
    SXFunction fcn("tmp", std::vector<SX>(), make_vector(ex));
    fcn.evaluate();
    return fcn.output();
  }

} // namespace casadi
