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

#include "multiplication.hpp"
#include "../matrix/matrix_tools.hpp"
#include <vector>

using namespace std;

namespace CasADi{

Multiplication::Multiplication(const MX& x, const MX& y){
  setDependencies(x,y);
  if(x.size2() != y.size1()) throw CasadiException("Multiplication::dimension mismatch");
  setSize(x.size1(),y.size2());
}

Multiplication* Multiplication::clone() const{
  return new Multiplication(*this);
}

void Multiplication::print(std::ostream &stream) const{
  stream << "prod(" << dep(0) << "," << dep(1) << ")";
}


void Multiplication::evaluate(int fsens_order, int asens_order){
  // Erase result
  fill(output().begin(), output().end(), 0);
  
  // Matrix multiplication
  matrix_matrix_mult(input(0),input(1),output());
  
  if(fsens_order>0){
    // Erase result
    fill(fwdSens().begin(), fwdSens().end(), 0);

    // Matrix multiplication, first argument
    matrix_matrix_mult(fwdSeed(0),input(1),fwdSens());

    // Matrix multiplication, second argument
    matrix_matrix_mult(input(0),fwdSeed(1),fwdSens());
  }
  
  if(asens_order>0){
    // Matrix multiplication, first argument
    matrix_matrix_trans_mult(adjSeed(),input(1),adjSens(0));

    // Matrix multiplication, second argument
    matrix_trans_matrix_mult(input(0),adjSeed(),adjSens(1));
  }
}


} // namespace CasADi

