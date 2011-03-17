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
typedef SX (*sder_ptr)(const SX&, const SX&, const SX&); // symbolic partial derivative
//@}
	
/** \brief  Zero derivative */
static SX zero_sder(const SX &f, const SX& x, const SX& y){ return 0; }

//@{
/** \brief  Addition */
static SX add_sder1(const SX &f, const SX& x, const SX& y){ return 1; }
static SX add_sder2(const SX &f, const SX& x, const SX& y){ return 1; }
//@}

//@{
/** \brief  Subtraction */
static SX sub_sder1(const SX &f, const SX& x, const SX& y){ return 1; }
static SX sub_sder2(const SX &f, const SX& x, const SX& y){ return -1; }
//@}

//@{
/** \brief  Multiplication */
static SX mul_sder1(const SX &f, const SX& x, const SX& y){ return y; }
static SX mul_sder2(const SX &f, const SX& x, const SX& y){ return x; }
//@}

//@{
/** \brief  Negation */
static SX neg_sder1(const SX &f, const SX& x, const SX& y){ return -1; }
//@}

//@{
/** \brief  Division */
static SX div_sder1(const SX &f, const SX& x, const SX& y){ return 1/y; }
static SX div_sder2(const SX &f, const SX& x, const SX& y){ return -f/y; }
//@}

//@{
/** \brief  Natural exponential */
static SX exp_sder1(const SX &f, const SX& x, const SX& y){ return f; }
//@}

//@{
/** \brief  Natural logarithm */
static SX log_sder1(const SX &f, const SX& x, const SX& y){ return 1/x; }
//@}

//@{
/** \brief  Power */
static SX pow_sder1(const SX &f, const SX& x, const SX& y){ return y*f/x; }
static SX pow_sder2(const SX &f, const SX& x, const SX& y){ return log(x)*f; }
//@}

//@{
/** \brief  Square root */
static SX sqrt_sder1(const SX &f, const SX& x, const SX& y){ return 1/(2*f);}
//@}

//@{
/** \brief  Sine */
static SX sin_sder1(const SX &f, const SX& x, const SX& y){ return cos(x); }
//@}

//@{
/** \brief  Cosine */
static SX cos_sder1(const SX &f, const SX& x, const SX& y){ return -sin(x); }
//@}

//@{
/** \brief  Tangens */
static SX tan_sder1(const SX &f, const SX& x, const SX& y){ SX cosx = cos(x); return 1/(cosx*cosx);}
//@}

//@{
/** \brief  arcsin */
static SX asin_sder1(const SX &f, const SX& x, const SX& y){ return 1/sqrt(1-x*x); }
//@}

//@{
/** \brief  arccos */
static SX acos_sder1(const SX &f, const SX& x, const SX& y){ return -1/sqrt(1-x*x);}
//@}

//@{
/** \brief  arctan */
static SX atan_sder1(const SX &f, const SX& x, const SX& y){ return 1/(1+x*x); }
//@}

//@{
/** \brief  Step node */
//@}

//@{
/** \brief  Floor */
//@}

//@{
/** \brief  Ceil */
//@}

//@{
/** \brief  Min */
static SX fmin_sder1(const SX &f, const SX& x, const SX& y){ return x<=y; }
static SX fmin_sder2(const SX &f, const SX& x, const SX& y){ return y<=x; }
//@}

//@{
/** \brief  Max */
static SX fmax_sder1(const SX &f, const SX& x, const SX& y){ return x>=y; }
static SX fmax_sder2(const SX &f, const SX& x, const SX& y){ return y>=x; }
//@}

//@{
/** \brief  Equality */
//@}

//@{
/** \brief  Error function */
static SX erf_sder1(const SX &f, const SX& x, const SX& y){ return (2/sqrt(M_PI))*exp(-x*x);}
//@}

///Array of symbolic partial derivatives
static sder_ptr sder1[] = { // 
  add_sder1, sub_sder1, mul_sder1, div_sder1, neg_sder1, exp_sder1, log_sder1, pow_sder1, sqrt_sder1, sin_sder1, cos_sder1, tan_sder1, asin_sder1, acos_sder1, atan_sder1, zero_sder,  zero_sder,   zero_sder,  zero_sder,      erf_sder1, fmin_sder1, fmax_sder1
};
/// Array of symbolic partial derivatives
static sder_ptr sder2[] = { 
  add_sder2, sub_sder2, mul_sder2, div_sder2, zero_sder, zero_sder, zero_sder, pow_sder2, zero_sder,  zero_sder, zero_sder, zero_sder, zero_sder,  zero_sder,  zero_sder,  zero_sder,  zero_sder,   zero_sder,  zero_sder,      zero_sder, fmin_sder2, fmax_sder2
};


} // namespace CasADi

#undef RES
#undef F
#undef FX
#undef FY

#endif //BINARY_FUNCTIONS_HPP
