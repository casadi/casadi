/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


/** \brief Demonstration of how to work with the error handling in CasADi
 * NOTE: Example is mainly intended for developers of CasADi.
 * CasADi provides a set of macros facilitating debugging. They are designed to
 * work in a similar way as the macros in "assert.h" in the C standard library
 * with the difference that the error message will be contained in a C++
 * exception rather than written to standard error and causing program termination.
 *
 * \author Joel Andersson
 * \date 2012
 */

#include "casadi/casadi.hpp"

bool bad_test(){
  try {
    // This will fail
    casadi_assert(false, "Failing assert");

    // Returns true, but the code won't reach this place
    return true;
  } catch (std::exception& e) {
    casadi_error("bad_test3 failed:\n" + std::string(e.what()));
  }
}

int main(){
  // Warning
  casadi_warning("This function will fail.");

  // Developer error
  casadi_assert(bad_test(), "");

  return 0;
}
