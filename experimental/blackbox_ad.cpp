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

// Generate an SX expression for a particular control flow
class Tracer : public GenericExpression<Tracer, bool> {
 public:
  // Nondifferentiated value
  double v;
  // SX expression
  SXElem s;
  // Constructor
  Tracer(double v, SXElem s) : v(v), s(s) {}
  // Constructor with implicit type conversion
  Tracer(double v) : v(v), s(v) {}
  // Binary operation
  static Tracer binary(casadi_int op, const Tracer& x, const Tracer& y) {
    // Perform operation numerically
    double f;
    casadi_math<double>::fun(op, x.v, y.v, f);
    // Perform operation symbolically
    SXElem f_sx;
    casadi_math<SXElem>::fun(op, x.s, y.s, f_sx);
    // Combine
    return Tracer(f, f_sx);
  }
  // Unary operation
  static Tracer unary(casadi_int op, const Tracer& x) {
    return binary(op, x, 0);
  }
  // Binary operation, boolean return
  static bool logic_binary(casadi_int op, const Tracer& x, const Tracer& y) {
    double f;
    casadi_math<double>::fun(op, x.v, y.v, f);
    return static_cast<bool>(f);  // type cast can be avoided by adding second template parameter
  }
  // Unary operation, boolean return
  static bool logic_unary(casadi_int op, const Tracer& x) {
    return logic_binary(op, x, 0);
  }
};

}  // namespace casadi

// Templated blackbox function
template<typename T>
T testfun(T x, T y) {
  // Conditionals
  T z = x > y ? x + sin(y) : y + sqrt(x);
  // While loop
  while (z < 20) {
    z = z + 4;
  }
  // Done
  return z;
}

void test_fwdad(double x0, double seed_x) {
  // Create an instance that contains the value of x and the forward seed
  casadi::FwdAD<1> x(x0, &seed_x);
  // Call templated function, calculating forward mode AD on the fly
  casadi::FwdAD<1> z = testfun(x, x);
  // Output result
  std::cout << "Operator overloading AD: value = " << z.v << ", sensitivity = "
    << z.s[0] << std::endl;
}

void test_tracer(double x0, double seed_x) {
  // Expression corresponding to CasADi function inputs
  casadi::SX x = casadi::SX::sym("x");
  // Get the nonzeros (corresponding to doubles)
  std::vector<casadi::SXElem> x_nz = x.nonzeros();
  // Create an instance that contains the value of x and the forward seed
  casadi::Tracer tx(x0, x_nz.at(0));
  // Call templated function, generating an SX expression along with the numerical value
  casadi::Tracer tz = testfun(tx, tx);
  double z0 = tz.v;
  casadi::SX z = tz.s;
  // Output result
  std::cout << "Tracing AD: value = " << z0 << ", expression = " << z << std::endl;
  // Generate SX Function
  casadi::Function f("f", {x}, {z}, {"x"}, {"z"});
  // Generate forward mode AD (or any other derivative)
  casadi::Function fwd_f = f.forward(1);
  // Evaluate derivative numerically
  double fwd_z;
  fwd_f({&x0, &z0, &seed_x}, {&fwd_z});
  // Output result
  std::cout << "Traced SX + SCT AD: sensitivity = " << fwd_z << std::endl;
}

int main(){
  // Perform memoryless operator overloading orward mode AD during evaluation
  test_fwdad(2, 1);
  // Generate an SX expression for the control flow, then use standard CasADi AD
  test_tracer(2, 1);

  return 0;
}
