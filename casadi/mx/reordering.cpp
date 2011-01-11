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

#include "reordering.hpp"
#include "../stl_vector_tools.hpp"

using namespace std;

namespace CasADi{

Reordering::Reordering(const MX &dep__){
  setDependencies(dep__);
}

void Reordering::evaluate(int fsens_order, int asens_order){
  
  if(fsens_order==0){
    for (int k=0;k<size1()*size2();k++) {
      output()[k]=input(k2l(k))[k2k(k)];
    }
  } else {
    for (int k=0;k<size1()*size2();k++) {
      fwdSens()[k]=fwdSeed(k2l(k))[k2k(k)];
    }
  }
  
  if(asens_order>0){
    for (int k=0;k<size1()*size2();k++) {
      adjSens(k2k(k))[k2l(k)]+=adjSeed()[k];
    }
  }
}

} // namespace CasADi
