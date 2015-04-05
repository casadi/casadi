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


#include "casadi_map.hpp"
#include "../function/function_internal.hpp"
#include "../std_vector_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#ifdef WITH_OPENMP
#include <omp.h>
#endif //WITH_OPENMP

using namespace std;

namespace casadi {

  Map::Map(const Function& fcn, const vector<vector<MX> >& arg) : fcn_(fcn) {
    // Number of calls
    n_ = arg.size();

    // Get all inputs
    int f_num_in = fcn_.getNumInputs();
    vector<MX> all_arg;
    all_arg.reserve(n_ * f_num_in);
    for (vector<vector<MX> >::const_iterator j=arg.begin(); j!=arg.end(); ++j) {
      casadi_assert(j->size()==n_);
      for (int i=0; i<n_; ++i) {
        if (j->at(i).shape()==fcn_.input(i).shape()) {
          // Insert sparsity projection nodes if needed
          all_arg.push_back(j->at(i).setSparse(fcn_.input(i).sparsity()));
        } else {
          // Different dimensions
          if (j->at(i).isEmpty() || fcn_.input(i).isEmpty()) { // NOTE: To permissive?
            // Replace nulls with zeros of the right dimension
            all_arg.push_back(MX::zeros(fcn_.input(i).sparsity()));
          } else if (j->at(i).isScalar()) {
            // Scalar argument means set all
            all_arg.push_back(MX(fcn_.input(i).sparsity(), j->at(i)));
          } else {
            // Mismatching dimensions
            casadi_error("Cannot create map node: Dimension mismatch for argument "
                         << i << ". Argument has shape " << j->at(i).shape()
                         << " but function input is " << fcn_.input(i).shape());
          }
        }
      }
    }
    casadi_assert(all_arg.size() == n_ * f_num_in);
    setDependencies(all_arg);
    setSparsity(Sparsity::scalar());
  }

  OmpMap::OmpMap(const Function& fcn, const vector<vector<MX> >& arg) : Map(fcn, arg) {
  }

  Map* Map::clone() const {
    return new Map(*this);
  }

  OmpMap* OmpMap::clone() const {
    return new OmpMap(*this);
  }

  void Map::printPart(ostream &stream, int part) const {
    if (part==0) {
      stream << "(";
    } else if (part==ndep()) {
      stream << ")";
    } else {
      stream << ",";
    }
  }

  void Map::evalD(cp_double* arg, p_double* res, int* itmp, double* rtmp) {
    int f_num_in = fcn_.getNumInputs();
    int f_num_out = fcn_.getNumOutputs();
    for (int i=0; i<n_; ++i) {
      fcn_->evalD(arg, res, itmp, rtmp);
      arg += f_num_in;
      res += f_num_out;
    }
  }

  void OmpMap::evalD(cp_double* arg, p_double* res, int* itmp, double* rtmp) {
#ifndef WITH_OPENMP
    // Not available, switching to serial mode
    Map::evalD(arg, res, itmp, rtmp);
#else // WITH_OPENMP
    int f_num_in = fcn_.getNumInputs();
    int f_num_out = fcn_.getNumOutputs();
    size_t ni, nr;
    fcn_.nTmp(ni, nr);
#pragma omp parallel for
    for (int i=0; i<n_; ++i) {
      fcn_->evalD(arg + i*f_num_in, res + i*f_num_out,
                  itmp + i*ni, rtmp + i*nr);
    }
#endif // WITH_OPENMP
  }

  int Map::nout() const {
    int f_num_out = fcn_.getNumOutputs();
    return n_ * f_num_out;
  }

  const Sparsity& Map::sparsity(int oind) const {
    int f_num_out = fcn_.getNumOutputs();
    return fcn_.output(oind % f_num_out).sparsity();
  }

  Function& Map::getFunction() {
    return fcn_;
  }

  void Map::evalSX(cp_SXElement* arg, p_SXElement* res, int* itmp, SXElement* rtmp) {
    int f_num_in = fcn_.getNumInputs();
    int f_num_out = fcn_.getNumOutputs();
    for (int i=0; i<n_; ++i) {
      fcn_->evalSX(arg, res, itmp, rtmp);
      arg += f_num_in;
      res += f_num_out;
    }
  }

  void Map::deepCopyMembers(map<SharedObjectNode*, SharedObject>& already_copied) {
    MXNode::deepCopyMembers(already_copied);
    fcn_ = deepcopy(fcn_, already_copied);
  }

  void Map::nTmp(size_t& ni, size_t& nr) {
    fcn_.nTmp(ni, nr);
  }

  void OmpMap::nTmp(size_t& ni, size_t& nr) {
    fcn_.nTmp(ni, nr);
    ni *= static_cast<size_t>(n_);
    nr *= static_cast<size_t>(n_);
  }

} // namespace casadi
