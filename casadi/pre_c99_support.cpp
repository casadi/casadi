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

#include "pre_c99_support.hpp"
#include <limits>
#include <algorithm>

#if __STDC_VERSION__ < 199901L
int isnan(double x) throw(){return x!=x;}
int isinf(double x) throw(){return isnan(x-x);}
double erf(double x) throw(){return std::numeric_limits<double>::quiet_NaN();}
double fmin(double x, double y) throw(){ return std::min(x,y);}
double fmax(double x, double y) throw(){ return std::max(x,y);}
#endif // __STDC_VERSION__ < 199901L

#ifndef WITHOUT_INT_MATH
double exp(int x) throw(){ return exp(double(x)); }
double log(int x) throw(){ return log(double(x)); }
double sqrt(int x) throw(){ return sqrt(double(x)); }
double pow(int x, int y) throw(){ return pow(double(x),double(y)); }
double sin(int x) throw(){ return sin(double(x)); }
double cos(int x) throw(){ return cos(double(x)); }
double tan(int x) throw(){ return tan(double(x)); }
double asin(int x) throw(){ return asin(double(x)); }
double acos(int x) throw(){ return acos(double(x)); }
double atan(int x) throw(){ return atan(double(x)); }
int floor(int x) throw(){ return x; }
int ceil(int x) throw(){ return x; }
#endif // WITHOUT_INT_MATH

namespace CasADi{
  double sign(double x) throw(){ 
    return x<0 ? -1 : x>0 ? 1 : x; // NOTE: sign(nan) == nan
  }
}


