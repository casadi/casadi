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

#ifndef BINARY_FUNCTIONS_HPP
#define BINARY_FUNCTIONS_HPP

#include "../elementary_functions.hpp"
#include "../pre_c99_support.hpp"

#define RES(i,j) res[i*(i+1)/2+j]
#define F RES(0,0)
#define FX RES(1,0)
#define FY RES(1,1)

namespace CasADi{
using namespace std;

//@{
/** \brief  Prototypes for a binary function */
typedef SX (*sfun_ptr)(const SX&, const SX&); // symbolic function
typedef SX (*sder_ptr)(const SX&, const SX&, const SX&); // symbolic partial derivative
typedef void (*nfun_ptr)(double, double, double*);
//@}
	
/** \brief  Zero derivative */
static SX zero_sder(const SX &f, const SX& x, const SX& y){ return 0; }

//@{
/** \brief  Addition */
static SX add_sfcn(const SX& x, const SX& y){ return x+y; }
static SX add_sder1(const SX &f, const SX& x, const SX& y){ return 1; }
static SX add_sder2(const SX &f, const SX& x, const SX& y){ return 1; }
static void add_nfcn0(double x, double y, double *res){ F = x+y;}
//@}

//@{
/** \brief  Subtraction */
static SX sub_sfcn(const SX& x, const SX& y){ return x-y; }
static SX sub_sder1(const SX &f, const SX& x, const SX& y){ return 1; }
static SX sub_sder2(const SX &f, const SX& x, const SX& y){ return -1; }
static void sub_nfcn0(double x, double y, double *res){ F = x-y;}
//@}

//@{
/** \brief  Multiplication */
static SX mul_sfcn(const SX& x, const SX& y){ return x*y; }
static SX mul_sder1(const SX &f, const SX& x, const SX& y){ return y; }
static SX mul_sder2(const SX &f, const SX& x, const SX& y){ return x; }
static void mul_nfcn0(double x, double y, double *res){ F = x*y;}
//@}

//@{
/** \brief  Negation */
static SX neg_sfcn(const SX& x, const SX& y){ return -x; }
static SX neg_sder1(const SX &f, const SX& x, const SX& y){ return -1; }
static void neg_nfcn0(double x, double y, double *res){ F = -x;}
//@}

//@{
/** \brief  Division */
static SX div_sfcn(const SX& x, const SX& y){ return x/y; }
static SX div_sder1(const SX &f, const SX& x, const SX& y){ return 1/y; }
static SX div_sder2(const SX &f, const SX& x, const SX& y){ return -f/y; }
static void div_nfcn0(double x, double y, double *res){ F = x/y;}
//@}

//@{
/** \brief  Natural exponential */
static SX exp_sfcn(const SX& x, const SX& y){ return exp(x); }
static SX exp_sder1(const SX &f, const SX& x, const SX& y){ return f; }
static void exp_nfcn0(double x, double y, double *res){ F = exp(x);}
//@}

//@{
/** \brief  Natural logarithm */
static SX log_sfcn(const SX& x, const SX& y){ return log(x); }
static SX log_sder1(const SX &f, const SX& x, const SX& y){ return 1/x; }
static void log_nfcn0(double x, double y, double *res){ F = log(x);}
//@}

//@{
/** \brief  Power */
static SX pow_sfcn(const SX& x, const SX& y){ return pow(x,y); }
static SX pow_sder1(const SX &f, const SX& x, const SX& y){ return y*f/x; }
static SX pow_sder2(const SX &f, const SX& x, const SX& y){ return log(x)*f; }
static void pow_nfcn0(double x, double y, double *res){ F = pow(x,y);}
//@}

//@{
/** \brief  Square root */
static SX sqrt_sfcn(const SX& x, const SX& y){ return sqrt(x); }
static SX sqrt_sder1(const SX &f, const SX& x, const SX& y){ return 1/(2*f);}
static void sqrt_nfcn0(double x, double y, double *res){ F = sqrt(x);}
//@}

//@{
/** \brief  Sine */
static SX sin_sfcn(const SX& x, const SX& y){ return sin(x); }
static SX sin_sder1(const SX &f, const SX& x, const SX& y){ return cos(x); }
static void sin_nfcn0(double x, double y, double *res){ F = sin(x);}
//@}

//@{
/** \brief  Cosine */
static SX cos_sfcn(const SX& x, const SX& y){ return cos(x); }
static SX cos_sder1(const SX &f, const SX& x, const SX& y){ return -sin(x); }
static void cos_nfcn0(double x, double y, double *res){ F = cos(x);}
//@}

//@{
/** \brief  Tangens */
static SX tan_sfcn(const SX& x, const SX& y){ return tan(x); }
static SX tan_sder1(const SX &f, const SX& x, const SX& y){ SX cosx = cos(x); return 1/(cosx*cosx);}
static void tan_nfcn0(double x, double y, double *res){ F = tan(x);}
//@}

//@{
/** \brief  arcsin */
static SX asin_sfcn(const SX& x, const SX& y){ return asin(x); }
static SX asin_sder1(const SX &f, const SX& x, const SX& y){ return 1/sqrt(1-x*x); }
static void asin_nfcn0(double x, double y, double *res){ F = asin(x);}
//@}

//@{
/** \brief  arccos */
static SX acos_sfcn(const SX& x, const SX& y){ return acos(x); }
static SX acos_sder1(const SX &f, const SX& x, const SX& y){ return -1/sqrt(1-x*x);}
static void acos_nfcn0(double x, double y, double *res){ F = acos(x);}
//@}

//@{
/** \brief  arctan */
static SX atan_sfcn(const SX& x, const SX& y){ return atan(x); }
static SX atan_sder1(const SX &f, const SX& x, const SX& y){ return 1/(1+x*x); }
static void atan_nfcn0(double x, double y, double *res){ F = atan(x);}
//@}

//@{
/** \brief  Step node */
static SX step_sfcn(const SX& x, const SX& y){ return x >= 0; }
static void step_nfcn0(double x, double y, double *res){ F = x >= 0;}
//@}

//@{
/** \brief  Floor */
static SX floor_sfcn(const SX& x, const SX& y){ return floor(x); }
static void floor_nfcn0(double x, double y, double *res){ F = floor(x);}
//@}

//@{
/** \brief  Ceil */
static SX ceil_sfcn(const SX& x, const SX& y){ return ceil(x); }
static void ceil_nfcn0(double x, double y, double *res){ F = ceil(x);}
//@}

//@{
/** \brief  Min */
static SX fmin_sfcn(const SX& x, const SX& y){ return fmin(x,y); }
static SX fmin_sder1(const SX &f, const SX& x, const SX& y){ return x<=y; }
static SX fmin_sder2(const SX &f, const SX& x, const SX& y){ return y<=x; }
static void fmin_nfcn0(double x, double y, double *res){ F = min(x,y);}
//@}

//@{
/** \brief  Max */
static SX fmax_sfcn(const SX& x, const SX& y){ return fmax(x,y); }
static SX fmax_sder1(const SX &f, const SX& x, const SX& y){ return x>=y; }
static SX fmax_sder2(const SX &f, const SX& x, const SX& y){ return y>=x; }
static void fmax_nfcn0(double x, double y, double *res){ F = max(x,y);}
//@}

//@{
/** \brief  Equality */
static SX equality_sfcn(const SX& x, const SX& y){ return x == y; }
static void equality_nfcn0(double x, double y, double *res){ F = x==y;}
//@}

//@{
/** \brief  Error function */
static SX erf_sfcn(const SX& x, const SX& y){ return erf(x); }
static SX erf_sder1(const SX &f, const SX& x, const SX& y){ return (2/sqrt(M_PI))*exp(-x*x);}
static void erf_nfcn0(double x, double y, double *res){ F = ERF(x);}
//@}

/// Vector of symbolic functions
static sfun_ptr sfcn[] = { // 
  add_sfcn,  sub_sfcn,  mul_sfcn,  div_sfcn,  neg_sfcn,  exp_sfcn,  log_sfcn,  pow_sfcn,  sqrt_sfcn,  sin_sfcn,  cos_sfcn,  tan_sfcn,  asin_sfcn,  acos_sfcn,  atan_sfcn,  step_sfcn,  floor_sfcn,  ceil_sfcn,  equality_sfcn,  erf_sfcn,  fmin_sfcn,  fmax_sfcn
};
///Array of symbolic partial derivatives
static sder_ptr sder1[] = { // 
  add_sder1, sub_sder1, mul_sder1, div_sder1, neg_sder1, exp_sder1, log_sder1, pow_sder1, sqrt_sder1, sin_sder1, cos_sder1, tan_sder1, asin_sder1, acos_sder1, atan_sder1, zero_sder,  zero_sder,   zero_sder,  zero_sder,      erf_sder1, fmin_sder1, fmax_sder1
};
/// Array of symbolic partial derivatives
static sder_ptr sder2[] = { 
  add_sder2, sub_sder2, mul_sder2, div_sder2, zero_sder, zero_sder, zero_sder, pow_sder2, zero_sder,  zero_sder, zero_sder, zero_sder, zero_sder,  zero_sder,  zero_sder,  zero_sder,  zero_sder,   zero_sder,  zero_sder,      zero_sder, fmin_sder2, fmax_sder2
};
/// Partial derivatives up to order 0
static nfun_ptr nfun0[] = {   
  add_nfcn0, sub_nfcn0, mul_nfcn0, div_nfcn0, neg_nfcn0, exp_nfcn0, log_nfcn0, pow_nfcn0, sqrt_nfcn0, sin_nfcn0, cos_nfcn0, tan_nfcn0, asin_nfcn0, acos_nfcn0, atan_nfcn0, step_nfcn0, floor_nfcn0, ceil_nfcn0, equality_nfcn0, erf_nfcn0, fmin_nfcn0, fmax_nfcn0
};


} // namespace CasADi

#undef RES
#undef F
#undef FX
#undef FY

#endif //BINARY_FUNCTIONS_HPP
