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
    int nret = (stop-start)/step + ((stop-start)%step>0);
    std::vector<int> ret(nret);
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
  
  Matrix<double> operator==(const Matrix<double>& a, const Matrix<double>& b){
    return a.binary_old(CasADi::casadi_operators<double>::equality, b);
  }
  
  Matrix<double> operator>=(const Matrix<double>& a, const Matrix<double>& b){
    Matrix<double> ret = a + b;
    if (a.sparsity()==b.sparsity()) {
      for (int i=0;i<a.size();++i) ret.at(i) = a.at(i) >= b.at(i);
    } else {
      casadi_assert_message(0,"not implemented");
    }
    return ret;
  }
  
  Matrix<double> operator<=(const Matrix<double>& a, const Matrix<double>& b){
    Matrix<double> ret = a + b;
    if (a.sparsity()==b.sparsity()) {
      for (int i=0;i<a.size();++i) ret.at(i) = a.at(i) <= b.at(i);
    } else {
      casadi_assert_message(0,"not implemented");
    }
    return ret;
  }

  
} // namespace CasADi

