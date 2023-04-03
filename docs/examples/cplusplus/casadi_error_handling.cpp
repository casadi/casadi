/*
 *    MIT No Attribution
 *
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
 *
 *    Permission is hereby granted, free of charge, to any person obtaining a copy of this
 *    software and associated documentation files (the "Software"), to deal in the Software
 *    without restriction, including without limitation the rights to use, copy, modify,
 *    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 *    permit persons to whom the Software is furnished to do so.
 *
 *    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 *    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 *    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
