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
template<int N>
class FwdAD : public GenericExpression<FwdAD<N>, bool> {
 public:
  // Nondifferentiated value
  double v;
  // Sensitivities
  double s[N];
  // Constructor with implicit type conversion
  FwdAD(double v) : v(v) {
    std::fill(s, s + N, 0);
  }
  // Constructor from arrays
  FwdAD(double v, double s[]) : v(v) {
    std::copy(s, s + N, this->s);
  }
  // Constructor with scalar multiplication of sensitivities
  FwdAD(double v, double a1, const double* s1) : v(v) {
    for (int k = 0; k < N; ++k) s[k] = a1 * s1[k];
  }
  // Constructor with linear combination of sensitivities
  FwdAD(double v, double a1, const double* s1, double a2, const double* s2) : v(v) {
    for (int k = 0; k < N; ++k) s[k] = a1 * s1[k] + a2 * s2[k];
  }
  // Binary operation
  static FwdAD binary(casadi_int op, const FwdAD& x, const FwdAD& y) {
    // Evaluate and get partial derivatives
    double d[2], f;
    casadi_math<double>::derF(op, x.v, y.v, f, d);
    // Return linear combination
    return FwdAD<N>(f, d[0], x.s, d[1], y.s);
  }
  // Unary operation
  static FwdAD unary(casadi_int op, const FwdAD& x) {
    // Evaluate and get partial derivatives
    double d[2], f;
    casadi_math<double>::derF(op, x.v, 0, f, d);
    // Return scalar multiplication
    return FwdAD<N>(f, d[0], x.s);
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
  casadi::FwdAD<1> x = 2;
  x.s[0] = 1;
  casadi::FwdAD<1> z = testfun(x, x);

  std::cout << "here: z.v = " << z.v << ", z.s = " << z.s[0] << std::endl;

  return 0;
}
