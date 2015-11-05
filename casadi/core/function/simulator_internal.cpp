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


#include "simulator_internal.hpp"
#include "ivpsol.hpp"
#include "../std_vector_tools.hpp"

using namespace std;
namespace casadi {


  SimulatorInternal::SimulatorInternal(const std::string& name, const Function& integrator)
    : FunctionInternal(name), integrator_(integrator) {
    grid_ = dynamic_cast<Ivpsol*>(integrator_.get())->grid_;
  }

  SimulatorInternal::~SimulatorInternal() {
  }

  void SimulatorInternal::init() {
    ischeme_ = Function::ivpsol_in();
    oscheme_ = Function::ivpsol_out();

    // Allocate inputs
    ibuf_.resize(get_n_in());
    for (int i=0; i<ibuf_.size(); ++i) {
      input(i) = DMatrix(get_sparsity_in(i));
    }

    // Allocate outputs
    obuf_.resize(get_n_out());
    for (int i=0; i<obuf_.size(); ++i) {
      output(i) = DMatrix(get_sparsity_out(i));
    }

    // Call base class method
    FunctionInternal::init();

    alloc(integrator_);
  }

  void SimulatorInternal::evalD(void* mem, const double** arg, double** res, int* iw, double* w) {
    // Arguments when calling the integrator
    double** res1 = res+n_out();
    copy_n(res, n_out(), res1);

    // Reset the integrator_
    dynamic_cast<Ivpsol*>(integrator_.get())
      ->reset(arg, res1, iw, w);

    // Advance solution in time
    for (int k=0; k<grid_.size(); ++k) {
      // Integrate to the output time
      dynamic_cast<Ivpsol*>(integrator_.get())->advance(k);

      // Save the outputs of the function
      for (int i=0; i<n_out(); ++i) {
        const Matrix<double> &res = integrator_.output(i);
        if (res1[i]) copy(res->begin(), res->end(), res1[i]);
        res1[i] += res.nnz();
      }
    }
  }

} // namespace casadi
