/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#include <iostream>
#include <fstream>
#include <string>

#include <kinetics.hpp>

using namespace std;

int main(int argc, char *argv[]){

try{
  
// SXMatrix t("t");
// 
//   SXMatrix alpha("alpha");
//   SXMatrix L("L");
//   
//   Frame f0("World Frame");
//   
//   SXMatrix a=TRz(alpha);
//   SXMatrix b=tr(0,0,L);
//   Frame f1("f1",f0,TRz(alpha)*tr(0,0,L));
// 
// 
//   // print the frame
//   cout << f1 << endl;
//   
  return 0;

} catch (const char * str){
  cerr << str << endl;
  return 1;
}
}
