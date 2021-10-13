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

// Class with forward derivative calculation
class FwdAD {
 public:
  // Nondifferentiated value
  double v;
  // Sensitivities
  double s;
  // Constructor with implicit type conversion
  FwdAD(double v, double s = 0) : v(v), s(s) {}
  // Addition
  inline friend FwdAD operator+(FwdAD x, FwdAD y) {
    return FwdAD(x.v + y.v, x.s + y.s);
  }
  // Sine
  inline friend FwdAD sin(FwdAD x) {
    return FwdAD(std::sin(x.v), std::cos(x.v) * x.s);
  }
  // Comparisons return booleans
  inline friend bool operator==(FwdAD x, FwdAD y) {
    return x.v == y.v;
  }
  inline friend bool operator!=(FwdAD x, FwdAD y) {
    return x.v != y.v;
  }
  inline friend bool operator<(FwdAD x, FwdAD y) {
    return x.v < y.v;
  }
  inline friend bool operator>(FwdAD x, FwdAD y) {
    return x.v > y.v;
  }
  inline friend bool operator<=(FwdAD x, FwdAD y) {
    return x.v <= y.v;
  }
  inline friend bool operator>=(FwdAD x, FwdAD y) {
    return x.v >= y.v;
  }
};


int main(){
  FwdAD x = 2;
  x.s = 1;
  FwdAD z = testfun(x, x);



  std::cout << "here: z.v = " << z.v << ", z.s = " << z.s << std::endl;

  return 0;
}
