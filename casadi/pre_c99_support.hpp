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

#ifndef PRE_C99_SUPPORT_HPP
#define PRE_C99_SUPPORT_HPP

// Get GCC version if GCC is used
#ifdef __GNUC__
#ifdef __GNUC_MINOR__
#ifdef __GNUC_PATCHLEVEL__
#define GCC_VERSION (__GNUC__ * 10000 +__GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif // __GNUC_PATCHLEVEL__
#endif // __GNUC_MINOR__
#endif // __GNUC__

// Disable some Visual studio warnings
#ifdef _MSC_VER

#pragma warning (disable:4996)

// warning C4018: '<' : signed/unsigned mismatch
#pragma warning (disable:4018)

// warning C4800: 'int' : forcing value to bool 'true'or 'false'(performance warning)
#pragma warning (disable:4800)
#endif

#if __STDC_VERSION__ < 199901L
// pre-C99
int isnan(double x) throw();
int isinf(double x) throw();
double erf(double x) throw();
double fmin(double x, double y) throw();
double fmax(double x, double y) throw();
#endif // __STDC_VERSION__ < 199901L

// Do we need to define elementary operations for int?
#ifdef GCC_VERSION
#if GCC_VERSION>=40200
#define WITHOUT_INT_MATH
#endif
#endif

// Resolve ambigous overloading
#ifndef WITHOUT_INT_MATH
namespace std{
  double exp(int x) throw();
  double log(int x) throw();
  double sqrt(int x) throw();
  double pow(int x, int y) throw();
  double sin(int x) throw();
  double cos(int x) throw();
  double tan(int x) throw();
  double asin(int x) throw();
  double acos(int x) throw();
  double atan(int x) throw();
  int floor(int x) throw();
  int ceil(int x) throw();
} // namespace std
#endif // WITHOUT_INT_MATH



#endif // PRE_C99_SUPPORT_HPP
