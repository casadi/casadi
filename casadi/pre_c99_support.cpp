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
double erf(double x) throw(){
  // Approximation found in sourceforge and modified, originally from numerical recepies in fortran
  double sx = x<0 ? -1 : x>0 ? 1 : x;
  double z = sx*x;
  double t = 1.0/(1.0+0.5*z);
  return 1.-sx*(t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
  t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
  t*(-0.82215223+t*0.17087277))))))))));
}
double fmin(double x, double y) throw(){ return std::min(x,y);}
double fmax(double x, double y) throw(){ return std::max(x,y);}
#endif // __STDC_VERSION__ < 199901L

namespace CasADi{
  double erfinv(double x) throw(){
    // Approximation found in sourceforge and modified: Not very efficent
    if(x>=1){
      return x==1 ? std::numeric_limits<double>::infinity() : std::numeric_limits<double>::quiet_NaN();
    } else if(x<=-1){
      return x==-1 ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::quiet_NaN();
    } else if(x<-0.7){
        double z = sqrt(-log((1.0+x)/2.0));
        return -(((1.641345311*z+3.429567803)*z-1.624906493)*z-1.970840454)/((1.637067800*z+3.543889200)*z+1.0);
    } else {
      double y;
      if(x<0.7){
        double z = x*x;
        y = x*(((-0.140543331*z+0.914624893)*z-1.645349621)*z+0.886226899)/((((-0.329097515*z+0.012229801)*z+1.442710462)*z-2.118377725)*z+1.0);
      } else {
        double z = sqrt(-log((1.0-x)/2.0));
        y = (((1.641345311*z+3.429567803)*z-1.624906493)*z-1.970840454)/((1.637067800*z+3.543889200)*z+1.0);
      }
      
      //polish x to full accuracy
      y = y - (erf(y) - x) / (2.0/sqrt(M_PI) * exp(-y*y));
      y = y - (erf(y) - x) / (2.0/sqrt(M_PI) * exp(-y*y));
      return y;
    }
  }
  
} // namespace CasADi


