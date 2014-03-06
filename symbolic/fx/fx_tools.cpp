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
         
    if (fabs(lb.at(i) - ub.at(i))<=tol) {
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
      indicator.at(0) = (fabs(v.at(i)-lb.at(i))<=tol)? 'X' : 'o';
      if (lb.at(i)==-std::numeric_limits<double>::infinity()) indicator.at(0)='8';

      indicator.at(indicator_length+1) = (fabs(v.at(i)-ub.at(i))<=tol)? 'X' : 'o';
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

FX parameterizeTime(FX dae) {

   // dimensionless time
   MX tau = MX::sym("tau");
   
   // The dae parameters augmented with t0 and tf
   MX P = MX::sym("P",1,2+dae.input(DAE_P).size());
   MX t0 = P[0];
   MX tf = P[1];
   
   std::vector<MX> dae_in(DAE_NUM_IN);
   std::vector<MX> dae_input = dae.symbolicInput();
   
   if (dae.input(DAE_T).size()==1) {
     dae_in[DAE_T]    = t0 + (tf-t0)*tau;
   }
   
   dae_in[DAE_P]    = reshape(trans(P[range(2,2+dae.input(DAE_P).size())]),dae.input(DAE_P).sparsity());
   dae_in[DAE_X]    = dae_input[DAE_X];

   std::vector<MX> ret_in(DAE_NUM_IN);
   ret_in[DAE_T]    = tau;
   ret_in[DAE_P]    = P;
   ret_in[DAE_X]    = dae_input[DAE_X];

   std::vector<MX> ret_out(DAE_NUM_OUT);
   ret_out[DAE_ODE] = (tf-t0)*dae.call(dae_in)[0];
   
   MXFunction ret(ret_in,ret_out);
   if (dae.isInit()) ret.init();
   
   // Expand if dae was an SXFunction
   if(is_a<SXFunction>(dae)){
     if(!ret.isInit()) ret.init();
     SXFunction ret_sx(ret);
     if(ret.isInit()) ret_sx.init();
     return ret_sx;
   } else {
     return ret;
   }
 }
 
FX parameterizeTimeOutput(FX f) {
  // dimensionless time
   MX tau = MX::sym("tau");
   
   // The f parameters augmented with t0 and tf
   MX P = MX::sym("P",1,2+f.input(DAE_P).size());
   MX t0 = P[0];
   MX tf = P[1];
   
   std::vector<MX> f_in(DAE_NUM_IN);
   std::vector<MX> f_input = f.symbolicInput();
   
   if (f.input(DAE_T).size()==1) {
     f_in[DAE_T]    = t0 + (tf-t0)*tau;
   }
   
   f_in[DAE_P]    = reshape(trans(P[range(2,2+f.input(DAE_P).size())]),f.input(DAE_P).sparsity());
   f_in[DAE_X]    = f_input[DAE_X];

   std::vector<MX> ret_in(DAE_NUM_IN);
   ret_in[DAE_T]    = tau;
   ret_in[DAE_P]    = P;
   ret_in[DAE_X]    = f_input[DAE_X];

   MXFunction ret(ret_in,f.call(f_in));
   
   if (f.isInit()) ret.init();
   
   // Expand if f was an SXFunction
   if(is_a<SXFunction>(f)){
     if(!ret.isInit()) ret.init();
     SXFunction ret_sx(ret);
     if(ret.isInit()) ret_sx.init();
     return ret_sx;
   } else {
     return ret;
   }
}

Matrix<double> numSample1D(FX &fx, const Matrix<double> &grid) {
  // Can be parallelized
  casadi_assert_message(fx.isInit(),"numSample1D:: supplied function must be initialized.");
  casadi_assert_message(fx.getNumInputs()>=1,"numSample1D:: supplied function must have at least one input.");
  casadi_assert_message(fx.input().size1()==1, "numSample1D:: supplied fx must have a row-matrix-like first input, but you supplied a shape " << fx.input().dimString() << ".");
  casadi_assert_message(fx.input().dense()==1, "numSample1D:: supplied fx must have dense input, but you supplied " << fx.input().dimString() << ".");
  casadi_assert_message(grid.size2()==fx.input().size2(), "numSample1D:: supplied grid has a shape " << grid.dimString() << ", but the row size does not match the row size of the supplied fx first input, which is " << fx.input().dimString() << ".");
  std::vector< Matrix<double> > ret(grid.size1());
  for (int j=0;j<grid.size1();++j) {
    fx.input().set(grid(j,ALL));
    fx.evaluate();
    ret[j] = Matrix<double>(fx.output());
  }
  return vertcat(ret);
}
    
Matrix<double> numSample1DT(FX &fx, const Matrix<double> &grid) {
  casadi_error("Not implemented yet");
}

Matrix<double> numSample2D(FX &fx, const Matrix<double> &grid) {
  casadi_error("Not implemented yet");
}


} // namespace CasADi

