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
int pre_c99_isnan(double x);
int pre_c99_isinf(double x);
double pre_c99_erf(double x);

#if __STDC_VERSION__ >= 199901L
// C99
#define ISNAN(x) isnan(x)
#define ISINF(x) isinf(x)
#define ERF(x) erf(x)
# else
// pre-C99
#define ISNAN(x) pre_c99_isnan(x)
#define ISINF(x) pre_c99_isinf(x)
#define ERF(x) pre_c99_erf(x)

#endif // __STDC_VERSION__ >= 199901L


#endif // PRE_C99_SUPPORT_HPP
