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


/** \brief Demonstration of sparsity propagation in CasADi
 * NOTE: Example is mainly intended for developers of CasADi.
 * This example demonstrates how dependencies can be propagated through symbolic expressions,
 * a low-level feature useful when determining sparsity patterns. This functionality is mostly
 * used internally in CasADi when calculating the sparsity pattern of Jacobians and Hessians.
 *
 * \author Joel Andersson
 * \date 2012
 */

#include "casadi/casadi.hpp"

using namespace casadi;

// Print a (typically 64-bit) unsigned integer
void print_binary(bvec_t v){
  for(int k=0; k<bvec_size; ++k){
    if(k%4==0) std::cout << " ";
    if(k%16==0) std::cout << " ";
    if(v & (bvec_t(1)<<k)){
      std::cout << 1;
    } else {
      std::cout << 0;
    }
  }
  std::cout << std::endl;
}

int main(){
  // Test both SX and MX
  for(int test=0; test<2; ++test){

    // Create a simple function
    Function f;
    if(test==0){
      std::cout << "SX:" << std::endl;
      SX x = SX::sym("x",3);
      SX z = x(0)*x(0)+x(2) + 3;
      f = Function("f", {x}, {z});
    } else {
      std::cout << "MX:" << std::endl;
      MX x = MX::sym("x",3);
      MX z = x(0)*x(0)+x(2) + 3;
      f = Function("f", {x}, {z});
    }

    // Memory for inputs and outputs
    std::vector<bvec_t> f_in(f.nnz_in(0), 0);
    std::vector<bvec_t> f_out(f.nnz_out(0), 0);

    // Propagate from input to output (forward mode)
    std::cout << "forward mode" << std::endl;

    // Make sure that the class is able to support the dependency propagation
    casadi_assert(f.has_spfwd(), "Forward sparsity propagation not supported");

    // Pass seeds
    f_in[0] = bvec_t(1) << 0; // seed in direction 0
    f_in[1] = bvec_t(1) << 2; // seed in direction 2
    f_in[2] = (bvec_t(1) << 4) | (bvec_t(1) << 63); // seed in direction 4 and 63

    // Reset sensitivities
    f_out[0] = 0;

    // Propagate dependencies
    f({&f_in.front()}, {&f_out.front()});

    // Print the result
    print_binary(f_out[0]);

    // Propagate from output to input (adjoint/reverse/backward mode)
    std::cout << "backward mode" << std::endl;

    // Make sure that the class is able to support the dependency propagation
    casadi_assert(f.has_sprev(), "Backward sparsity propagation not supported");

    // Pass seeds
    f_out[0] = (bvec_t(1) << 5) | (bvec_t(1) << 6); // seed in direction 5 and 6

    // Reset sensitivities
    f_in[0] = 0;
    f_in[1] = 0;
    f_in[2] = 0;

    // Propagate dependencies
    f.rev({&f_in.front()}, {&f_out.front()});

    // Print the result
    print_binary(f_in[0]);
    print_binary(f_in[1]);
    print_binary(f_in[2]);
  }

  return 0;
}
