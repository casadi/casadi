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

#include "casadi/casadi.hpp"

using namespace std;
using namespace CasADi;

int main() {
  // Create an object that can be referenced weakly (deriving from CachedObject)
  CRSSparsity a = sp_diag(3);
  cout << "a = " << a << endl;
  
  // Create a weak reference
  WeakRef a_weak = a;
  
  // Check if the weak reference is null
  cout << "a_weak.isNull() = " << a_weak.isNull() << " ( == 0 since a is an owner)" << endl;
  
  // Get a copy of the object from the weak reference
  CRSSparsity a_copy = a_weak;
  cout << "a_copy = " << a_copy << endl;
  
  // Remove the owner
  a = sp_diag(4);
  
  // Check if the weak reference has become null
  cout << "a_weak.isNull() = " << a_weak.isNull() << " ( == 0 since a_copy is an owner)" << endl;
  
  // Remove the copy
  a_copy = sp_diag(5);
  
  // Check if the weak reference has become null
  cout << "a_weak.isNull() = " << a_weak.isNull() << " ( == 1 since both a and a_copy has been overwritten, leaving no owner)" << endl;
  
  return 0;
}
