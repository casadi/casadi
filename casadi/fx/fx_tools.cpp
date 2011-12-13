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

#include "fx_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../stl_vector_tools.hpp"
#include "integrator.hpp"
#include <iostream>

using namespace std;

namespace CasADi{
    
void reportConstraints(std::ostream &stream,const Matrix<double> &v, const Matrix<double> &lb, const Matrix<double> &ub, const std::string &name, double tol) { 

  casadi_assert_message(v.sparsity()==lb.sparsity(),"reportConstraints:: sparsities must match");
  casadi_assert_message(ub.sparsity()==lb.sparsity(),"reportConstraints:: sparsities must match");
  
  // Backup the formatting
  ios::fmtflags fmtflags_backup = stream.flags();
  int streamsize_backup = stream.precision();
  
  //stream.setf(ios::fixed,ios::floatfield);
  //stream.precision(8);
  
  // Check if any constraint is violated
  if ( all(v <= ub + tol) && all(v >= lb - tol) ) {
    stream << "All " << v.size() << " constraints on " << name << " are met: " << endl;
  } else {
    stream << "Problem with constraints on " << name << ": " << endl;
  }
  
  // Make a horizontal rule
  stream.width(60);
  stream.fill('-');
  stream  << "-" << endl;
  stream.fill(' ');

  // The length of the numeric fields
  int fieldlength = 10;
  // The length of the constraint visualizer strip
  int indicator_length = 15;
      
  // Loop over the elements of v
  for (int i=0;i<v.size();i++) {
  
    stream.width(5);
    stream << i << ". |   ";
         
    if (abs(lb.at(i) - ub.at(i))<=tol) {
       stream.width(fieldlength);
       stream << lb.at(i) << " ==  ";
       stream.width(fieldlength);
       stream << v.at(i) << "      ";
       stream.width(fieldlength);
       stream << " ";
       stream << "  |   ";
    } else {
      // BEGIN  - construct the constraint visualizer strip
      std::string indicator(indicator_length+2,'-');
      indicator.at(0) = (abs(v.at(i)-lb.at(i))<=tol)? 'X' : 'o';
      if (lb.at(i)==-std::numeric_limits<double>::infinity()) indicator.at(0)='8';

      indicator.at(indicator_length+1) = (abs(v.at(i)-ub.at(i))<=tol)? 'X' : 'o';
      if (ub.at(i)==std::numeric_limits<double>::infinity()) indicator.at(indicator_length+1)='8';
            
      if (v.at(i) <= (ub.at(i) + tol) && v.at(i) >= (lb.at(i) - tol)) {
        int index = (v.at(i)-lb.at(i))/(ub.at(i)-lb.at(i))*(indicator_length-1);
        index = min(max(0,index),indicator_length-1);
        indicator.at(1+index) = '=';
      }
      // END - construct the constraint visualizer strip
      
      stream.width(fieldlength);
      stream << lb.at(i) << " <=  ";
      stream.width(fieldlength);
      stream << v.at(i) << " <= ";
      stream.width(fieldlength);
      stream << ub.at(i) << "    | ";
      if (v.at(i) <= (ub.at(i) + tol) && v.at(i) >= (lb.at(i) - tol)) {
        stream  << indicator;
      }
    }

    if (v.at(i) <= (ub.at(i) + tol) && v.at(i) >= (lb.at(i) - tol)) {
    } else {
        stream  << "  VIOLATED";
    }
    stream  << endl;
  }
  
  // Make a horizontal rule
  stream.width(60);
  stream.fill('-');
  stream  << "-" << endl;
  stream.fill(' ');
  
  // Restore the formatting
  stream.setf(fmtflags_backup);
  stream.precision(streamsize_backup);
}

MXFunction parameterizeTime(FX &dae) {
   casadi_warning("This piece of code lacks a unittest");
   
   // dimensionless time
   MX tau("tau");
   
   // The dae parameters augmented with t0 and tf
   MX P("P",2+dae.input(DAE_P).size());
   MX t0 = P[0];
   MX tf = P[1];
   
   std::vector<MX> dae_in(DAE_NUM_IN);
   if (dae.input(DAE_T).size()==1) {
     dae_in[DAE_T]    = t0 + (tf-t0)*tau;
   }
   dae_in[DAE_P]    = reshape(P[range(2,dae.input(DAE_P).size())],dae.input(DAE_P).sparsity());
   dae_in[DAE_Y]    = dae.input(DAE_Y);
   dae_in[DAE_YDOT] = dae.input(DAE_YDOT);

   std::vector<MX> ret_in(DAE_NUM_IN);
   ret_in[DAE_T]    = tau;
   ret_in[DAE_P]    = P;
   ret_in[DAE_Y]    = dae.input(DAE_Y);
   ret_in[DAE_YDOT] = dae.input(DAE_YDOT);

   std::vector<MX> ret_out(DAE_NUM_OUT);
   ret_out[DAE_RES] = (tf-t0)*dae.call(dae_in)[0];
   
   MXFunction ret(ret_in,ret_out);
   
   if (dae.isInit()) ret.init();
   
   return ret;
   
 }


} // namespace CasADi

