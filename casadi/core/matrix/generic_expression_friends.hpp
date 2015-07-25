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


/** \brief This file is included from generic_expression.hpp and used for SWIG wrapping
    \author Joel Andersson
    \date 2015
*/

/// Addition
inline SWIG_FRIEND ExType operator+(const ExType &x, const ExType &y) { return x.zz_plus(y); }

/// Subtraction
inline SWIG_FRIEND ExType operator-(const ExType &x, const ExType &y) { return x.zz_minus(y); }

/// Elementwise multiplication
inline SWIG_FRIEND ExType operator*(const ExType &x, const ExType &y) { return x.zz_times(y); }

/// Elementwise division
inline SWIG_FRIEND ExType operator/(const ExType &x, const ExType &y) { return x.zz_rdivide(y); }

/// Logic less than
inline SWIG_FRIEND ExType operator<(const ExType &x, const ExType &y) { return x.zz_lt(y); }

/// Logic less or equal to
inline SWIG_FRIEND ExType operator<=(const ExType &x, const ExType &y) { return x.zz_le(y); }

/// Logic greater than
inline SWIG_FRIEND ExType operator>(const ExType &x, const ExType &y) { return x.zz_gt(y); }

/// Logic greater or equal to
inline SWIG_FRIEND ExType operator>=(const ExType &x, const ExType &y) { return x.zz_ge(y); }

/// Logic equal to
inline SWIG_FRIEND ExType operator==(const ExType &x, const ExType &y) { return x.zz_eq(y); }

/// Logic not equal to
inline SWIG_FRIEND ExType operator!=(const ExType &x, const ExType &y) { return x.zz_ne(y); }

/// Logic and
inline SWIG_FRIEND ExType operator&&(const ExType &x, const ExType &y) { return x.zz_and(y); }

/// Logic or
inline SWIG_FRIEND ExType operator||(const ExType &x, const ExType &y) { return x.zz_or(y); }

/** \brief  Simplify an expression */
inline SWIG_FRIEND ExType simplify(const ExType &x) { return x.zz_simplify();}

/** \brief Check if two nodes are equivalent up to a given depth.
 *  Depth=0 checks if the expressions are identical, i.e. points to the same node.
 *
 *  a = x*x
 *  b = x*x
 *
 *  a.isEqual(b, 0)  will return false, but a.isEqual(b, 1) will return true
 */
inline SWIG_FRIEND bool isEqual(const ExType& x, const ExType& y, int depth=0) {
  return x.zz_isEqual(y, depth);
}

/** \brief  check if the matrix is 0 (note that false negative answers are possible) */
inline SWIG_FRIEND bool iszero(const ExType& x) { return x.isZero();}

/** \brief Absolute value, C++ syntax */
inline SWIG_FRIEND ExType abs(const ExType& x) { return x.zz_abs();}

/** \brief Absolute value, C syntax */
inline SWIG_FRIEND ExType fabs(const ExType& x) { return x.zz_abs();}

/** \brief Square root */
inline SWIG_FRIEND ExType sqrt(const ExType& x) { return x.zz_sqrt();}

/** \brief Sine */
inline SWIG_FRIEND ExType sin(const ExType& x) { return x.zz_sin();}

/** \brief Cosine */
inline SWIG_FRIEND ExType cos(const ExType& x) { return x.zz_cos();}

/** \brief Tangent */
inline SWIG_FRIEND ExType tan(const ExType& x) { return x.zz_tan();}

/** \brief Arc tangent */
inline SWIG_FRIEND ExType atan(const ExType& x) { return x.zz_atan();}

/** \brief Arc sine */
inline SWIG_FRIEND ExType asin(const ExType& x) { return x.zz_asin();}

/** \brief Arc cosine */
inline SWIG_FRIEND ExType acos(const ExType& x) { return x.zz_acos();}

/** \brief Hyperbolic tangent */
inline SWIG_FRIEND ExType tanh(const ExType& x) { return x.zz_tanh();}

/** \brief Hyperbolic sine */
inline SWIG_FRIEND ExType sinh(const ExType& x) { return x.zz_sinh();}

/** \brief Hyperbolic cosine */
inline SWIG_FRIEND ExType cosh(const ExType& x) { return x.zz_cosh();}

/** \brief Arc hyperbolic tangent */
inline SWIG_FRIEND ExType atanh(const ExType& x) { return x.zz_atanh();}

/** \brief Arc hyperbolic sine */
inline SWIG_FRIEND ExType asinh(const ExType& x) { return x.zz_asinh();}

/** \brief Arc hyperbolic cosine */
inline SWIG_FRIEND ExType acosh(const ExType& x) { return x.zz_acosh();}

/** \brief Natural exponential function (elementwise for matrix types) */
inline SWIG_FRIEND ExType exp(const ExType& x) { return x.zz_exp();}

/** \brief Natural logarithm */
inline SWIG_FRIEND ExType log(const ExType& x) { return x.zz_log();}

/** \brief 10-base logarithm */
inline SWIG_FRIEND ExType log10(const ExType& x) { return x.zz_log10();}

/** \brief Round down to nearest integer */
inline SWIG_FRIEND ExType floor(const ExType& x) { return x.zz_floor();}

/** \brief Round up to nearest integer */
inline SWIG_FRIEND ExType ceil(const ExType& x) { return x.zz_ceil();}

/** \brief Error function */
inline SWIG_FRIEND ExType erf(const ExType& x) { return x.zz_erf();}

/** \brief Sign function (note sign(nan) == nan, sign(0) == 0) */
inline SWIG_FRIEND ExType sign(const ExType& x) { return x.zz_sign();}

/** \brief Power (elementwise for matrix types) */
inline SWIG_FRIEND ExType pow(const ExType& x, const ExType& n) { return x.zz_power(n);}

/** \brief Modulo */
inline SWIG_FRIEND ExType fmod(const ExType& x, const ExType& y) { return x.zz_mod(y);}

/** \brief Arctan2 */
inline SWIG_FRIEND ExType atan2(const ExType& x, const ExType& y) { return x.zz_atan2(y);}

/** \brief Minimum of two values */
inline SWIG_FRIEND ExType fmin(const ExType& x, const ExType& y) { return x.zz_min(y);}

/** \brief Maximum of two values */
inline SWIG_FRIEND ExType fmax(const ExType& x, const ExType& y) { return x.zz_max(y);}
