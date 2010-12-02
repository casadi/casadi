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

#ifndef SUPERLU_INTERNAL_HPP
#define SUPERLU_INTERNAL_HPP

#include "superlu.hpp"
#include "external_packages/superlu_4_1/SRC/slu_ddefs.h"

namespace CasADi{
  
class SuperLUInternal : public FXNode{
  public:
    // Create a linear solver given a sparsity pattern and a number of right hand sides
    SuperLUInternal(int nrow, int ncol, const std::vector<int>& rowind, const std::vector<int>& col, int nrhs);

    // Destructor
    virtual ~SuperLUInternal();
    
    // Initialize the solver
    virtual void init();

    // Solve the system of equations
    virtual void evaluate(int fsens_order, int asens_order);

    // Solve the system of equations
    void solve(Factorization fact);

  protected:
    
    // Sparsity in CRS format
    int nrow_, ncol_;
    std::vector<int> rowind_, col_;
    int nrhs_;
    
    // Is initialized
    bool is_init;
    
    // SuperLU data structures
    SuperMatrix A, L, U, B;

    double   *a, *rhs;
    int      *asub, *xa;
    int      *perm_r; /* row permutations from partial pivoting */
    int      *perm_c; /* column permutation vector */
    int      info, i, permc_spec;
    superlu_options_t options;
    SuperLUStat_t stat;

    
};

} // namespace CasADi

#endif //SUPERLU_INTERNAL_HPP

