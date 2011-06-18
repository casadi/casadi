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

#include "matrix_tools.hpp"

namespace CasADi{
  
  std::vector<int> range(int start, int stop, int step, int len){
    start = std::min(start,len);
    stop = std::min(stop,len);
    std::vector<int> ret((stop-start)/step);
    int ind = start;
    for(std::vector<int>::iterator it=ret.begin(); it!=ret.end(); ++it){
      *it = ind;
      ind += step;
    }
    return ret;
  }
  
  std::vector<int> range(int stop){
    return range(0,stop);
  }
  
  
  Matrix<double> operator==(Matrix<double>& a, Matrix<double>& b) {
    std::cout << "Welcome, stranger" << std::endl;
    
    if (a.numel()==1) {
      double a_ = (a.size()==1) ? a.at(0) : 0;
      
      // If b is scalar, do a quick check
      if (b.numel()==1 && b.size()==1) return a_==b.at(0);
        
      // If b is sparse and a non-zero, then return false
      if (a_!=0 && b.numel()!=b.size()) return false;
      
      throw CasadiException("operator==: Not implemented");
      
      // Else, iterate over all elements b
      for (int k=0;k<b.size();k++) {
        if (b.data()[k]!=a_) return false;
      }
      return true;
    } else if (b.numel()==1) {
      throw CasadiException("operator==: Not implemented");
      /*return operator==(b,a);*/
    } else if (a.size1()==b.size1() && a.size2()==b.size2()) {
      throw CasadiException("operator==: Not implemented");
      /*std::vector<int> mapping;
      CRSSparsity sp = a.sparsity().patternUnion(a.sparsity(),mapping);
      int elA=0, elB=0;
      for(int k=0; k<mapping.size(); ++k){
        if(mapping[k]<0) {
           if (0!=a.at(elA++)) return false;
        } else if(mapping[k]>0){
           if (0!=b.at(elB++)) return false;
        } else {
          if (a.at(elB++)!=b.at(elB++)) return false;
        }
      }
      return true;*/
    } else {
      std::stringstream ss;
      ss << "operator==:Dimensions mismatch. a and b in a==b must be of same shape" << std::endl;
      ss << "Got a=" << a.dimString() <<  std::endl;
      ss << "Got b=" << a.dimString() <<  std::endl;
      throw CasadiException(ss.str());
    }
  }

} // namespace CasADi

