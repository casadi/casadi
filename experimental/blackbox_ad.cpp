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

#include <casadi/casadi.hpp>

#include <iomanip>
#include <ctime>
#include <cstdlib>

namespace casadi {

// Class with forward derivative calculation
class FwdAD : public GenericExpression<FwdAD, bool> {
 public:
  // Nondifferentiated value
  double v;
  // Sensitivities
  double s;
  // Constructor with implicit type conversion
  FwdAD(double v, double s = 0) : v(v), s(s) {}
  // Binary operation
  static FwdAD binary(casadi_int op, const FwdAD& x, const FwdAD& y) {
    // Evaluate and get partial derivatives
    double d[2], f;
    casadi_math<double>::derF(op, x.v, y.v, f, d);
    // Propagate
    return FwdAD(f, x.s * d[0] + y.s * d[1]);
  }
  // Unary operation
  static FwdAD unary(casadi_int op, const FwdAD& x) {
    return binary(op, x, 0);
  }
  // Binary operation, boolean return
  static bool logic_binary(casadi_int op, const FwdAD& x, const FwdAD& y) {
    double f;
    casadi_math<double>::fun(op, x.v, y.v, f);
    return static_cast<bool>(f);  // type cast can be avoided by adding second template parameter
  }
  // Unary operation, boolean return
  static bool logic_unary(casadi_int op, const FwdAD& x) {
    return logic_binary(op, x, 0);
  }
};

}  // namespace casadi

// Templated blackbox function
template<typename T>
T testfun(T x, T y) {
  // Conditionals
  T z = x > y ? x + sin(y) : y + sin(x);
  // While loop
  while (z < 20) {
    z = z + 4;
  }
  // Done
  return z;
}

int main(){
  casadi::FwdAD x = 2;
  x.s = 1;
  casadi::FwdAD z = testfun(x, x);

  std::cout << "here: z.v = " << z.v << ", z.s = " << z.s << std::endl;

  return 0;
}
