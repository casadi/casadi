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

  Map::Map(const Function& fcn, const vector< MX >& arg) :
      fcn_(fcn), n_in_(fcn.nIn()), n_out_(fcn.nOut()) {

    casadi_assert(fcn_.nIn()==arg.size());

    step_in_.resize(n_in_, 0);
    step_out_.resize(n_out_, 0);

    // Get all inputs
    int f_num_in = n_in_;
    vector<MX> all_arg;

    // Number of calls
    n_ = 1;
    for (int i=0;i<n_in_;++i) {
      casadi_assert(arg[i].size1()==fcn_.input(i).size1());
      casadi_assert(arg[i].size2() % fcn_.input(i).size2()==0);
      int n = arg[i].size2() / fcn_.input(i).size2();

      // An input is either of original shape or repeated n times
      casadi_assert(n==1 || n_==1 || n==n_);
      if (n>1) n_ = n;
      all_arg.push_back(project(arg[i], repmat(fcn_.input(i).sparsity(), 1, n)));

      if (n!=1) {
        step_in_[i] = fcn_.input(i).nnz();
      }
    }

    output_sparsity_.resize(n_out_);
    for (int i=0;i<n_out_;++i) {
      output_sparsity_[i] = repmat(fcn_.output(i).sparsity(), 1, n_);
      step_out_[i] = fcn_.output(i).nnz();
    }

    setDependencies(all_arg);
    setSparsity(Sparsity::scalar());
  }

  OmpMap::OmpMap(const Function& fcn, const vector<MX >& arg) : Map(fcn, arg) {
  }

  Map* Map::clone() const {
    return new Map(*this);
  }

  OmpMap* OmpMap::clone() const {
    return new OmpMap(*this);
  }

  std::string Map::print(const std::vector<std::string>& arg) const {
    stringstream ss;
    ss << fcn_.getOption("name") << ".map(";
    for (int i=0; i<ndep(); ++i) {
      if (i!=0) ss << ", ";
      ss << arg.at(i);
    }
    ss << ")";
    return ss.str();
  }

  void Map::evalD(const double** arg, double** res, int* iw, double* w) {
    const double** arg1 = arg+ndep();
    double** res1 = res+nout();

    for (int i=0; i<n_; ++i) {
      for (int j=0; j<n_in_; ++j) {
        arg1[j] = arg[j]+i*step_in_[j];
      }
      for (int j=0; j<n_out_; ++j) {
        res1[j] = (res[j]==0)? 0: res[j]+i*step_out_[j];
      }
      fcn_->evalD(arg1, res1, iw, w);
    }
  }

  void OmpMap::evalD(const double** arg, double** res, int* iw, double* w) {
#ifndef WITH_OPENMP
    // Not available, switching to serial mode
    Map::evalD(arg, res, iw, w);
#else // WITH_OPENMP
    size_t sz_arg, sz_res, sz_iw, sz_w;
    fcn_.sz_work(sz_arg, sz_res, sz_iw, sz_w);

#pragma omp parallel for
    for (int i=0; i<n_; ++i) {
      const double** arg_i = arg + n_in_ + sz_arg*i;
      for (int j=0; j<n_in_; ++j) {
        arg_i[j] = arg[j]+i*step_in_[j];
      }
      double** res_i = res + n_out_ + sz_res*i;
      for (int j=0; j<n_out_; ++j) {
        res_i[j] = (res[j]==0)? 0: res[j]+i*step_out_[j];
      }
      int* iw_i = iw + i*sz_iw;
      double* w_i = w + i*sz_w;
      fcn_->evalD(arg_i, res_i, iw_i, w_i);
    }
#endif // WITH_OPENMP
  }

  void Map::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    const bvec_t** arg1 = arg+ndep();
    bvec_t** res1 = res+nout();
    for (int i=0; i<n_; ++i) {
      for (int j=0; j<n_in_; ++j) {
        arg1[j] = arg[j]+i*step_in_[j];
      }
      for (int j=0; j<n_out_; ++j) {
        res1[j] = (res[j]==0)? 0: res[j]+i*step_out_[j];
      }
      fcn_->spFwd(arg1, res1, iw, w);
    }
  }

  void Map::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    bvec_t** arg1 = arg+ndep();
    bvec_t** res1 = res+nout();
    for (int i=0; i<n_; ++i) {
      for (int j=0; j<n_in_; ++j) {
        arg1[j] = arg[j]+i*step_in_[j];
      }
      for (int j=0; j<n_out_; ++j) {
        res1[j] = (res[j]==0)? 0: res[j]+i*step_out_[j];
      }
      fcn_->spAdj(arg1, res1, iw, w);
    }
  }

  int Map::nout() const {
    return n_out_;
  }

  const Sparsity& Map::sparsity(int oind) const {
    return output_sparsity_[oind];
  }

  const Function& Map::getFunction() const {
    return fcn_;
  }

  void Map::evalSX(const SXElement** arg, SXElement** res, int* iw, SXElement* w) {
    const SXElement** arg1 = arg+ndep();
    SXElement** res1 = res+nout();
    for (int i=0; i<n_; ++i) {
      for (int j=0; j<n_in_; ++j) {
        arg1[j] = arg[j]+i*step_in_[j];
      }
      for (int j=0; j<n_out_; ++j) {
        res1[j] = (res[j]==0)? 0: res[j]+i*step_out_[j];
      }
      fcn_->evalSX(arg1, res1, iw, w);
    }
  }

  void Map::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    // Collect arguments
    vector<vector<MX> > v(n_);
    vector<MX>::const_iterator arg_it = arg.begin();
    for (int i=0; i<n_; ++i) {
      v[i] = vector<MX>(arg_it, arg_it+n_in_);
      arg_it += n_in_;
    }

    // Call in parallel
    v = fcn_.map(v, parallelization());

    // Get results
    vector<MX>::iterator res_it = res.begin();
    for (int i=0; i<n_; ++i) {
      copy(v[i].begin(), v[i].end(), res_it);
      res_it += n_out_;
    }
  }

  void Map::evalFwd(const std::vector<std::vector<MX> >& fseed,
                     std::vector<std::vector<MX> >& fsens) {

    // Derivative function
    int nfwd = fsens.size();

    Function dfcn = fcn_.derForward(nfwd);

    // Nondifferentiated inputs and outputs
    vector<MX> arg(ndep());
    for (int i=0; i<arg.size(); ++i) arg[i] = dep(i);
    vector<MX> res(nout());
    for (int i=0; i<res.size(); ++i) res[i] = getOutput(i);

    // Collect arguments
    vector< MX > v;
    v.insert(v.end(), arg.begin(), arg.end());
    v.insert(v.end(), res.begin(), res.end());
    for (int j=0;j<nfwd;++j) {
      v.insert(v.end(), fseed[j].begin(), fseed[j].end());
    }

    // Call the cached function
    v = dfcn.map(v, parallelization());

    for (int k=0;k<n_out_;++k) {
      for (int j=0;j<nfwd;++j) {
        fsens[j][k] = v[k+j*n_out_];
      }
    }
  }

  void Map::evalAdj(const std::vector<std::vector<MX> >& aseed,
                     std::vector<std::vector<MX> >& asens) {
    // Derivative function
    int nadj = asens.size();
    Function dfcn = fcn_.derReverse(nadj);

    // Nondifferentiated inputs and outputs
    vector<MX> arg(ndep());
    for (int i=0; i<arg.size(); ++i) arg[i] = dep(i);
    vector<MX> res(nout());
    for (int i=0; i<res.size(); ++i) res[i] = getOutput(i);

    // Collect arguments
    vector< MX > v;
    v.insert(v.end(), arg.begin(), arg.end());
    v.insert(v.end(), res.begin(), res.end());
    for (int j=0;j<nadj;++j) {
      v.insert(v.end(), aseed[j].begin(), aseed[j].end());
    }

    // Call the cached function
    v = dfcn.map(v, parallelization());
    for (int k=0;k<n_in_;++k) {
      for (int j=0;j<nadj;++j) {
        if (v[k+j*n_in_].size2()>asens[j][k].size2()) {
          // Deal with non-repeated inputs
          vector< MX > s = horzsplit(v[k+j*n_in_], asens[j][k].size2());
          for (int i=0;i<s.size();++i) {
            asens[j][k] += s[i];
          }
        } else {
          asens[j][k] += v[k+j*n_in_];
        }
      }
    }
  }

  void Map::deepCopyMembers(map<SharedObjectNode*, SharedObject>& already_copied) {
    MXNode::deepCopyMembers(already_copied);
    fcn_ = deepcopy(fcn_, already_copied);
  }

  size_t Map::sz_arg() const {
    return ndep() + fcn_.sz_arg();
  }

  size_t Map::sz_res() const {
    return nout() + fcn_.sz_res();
  }

  size_t Map::sz_iw() const {
    return fcn_.sz_iw();
  }

  size_t Map::sz_w() const {
    return fcn_.sz_w();
  }

  size_t OmpMap::sz_arg() const {
    return ndep() + fcn_.sz_arg()*n_;
  }

  size_t OmpMap::sz_res() const {
    return nout() + fcn_.sz_res()*n_;
  }

  size_t OmpMap::sz_iw() const {
    return fcn_.sz_iw()*n_;
  }

  size_t OmpMap::sz_w() const {
    return fcn_.sz_w()*n_;
  }

  std::vector<MX >
  Map::create(const Function& fcn, const std::vector< MX > &arg,
              const std::string& parallelization) {

    if (parallelization.compare("expand")==0) {
      //Bypass the Map, call the original function n times
      int n = 1;

      // split arguments
      std::vector< std::vector< MX > > arg_split;
      for (int i=0;i<arg.size();++i) {
        arg_split.push_back(horzsplit(arg[i], fcn.input(i).size2()));
        if (arg_split[i].size()>n) n = arg_split[i].size();
      }

      std::vector< std::vector<MX> > ret(n);

      // call the original function
      for (int i=0; i<n; ++i) {
        std::vector< MX > argi;
        for (int j=0;j<arg.size();++j) {
          if (arg_split[j].size()>1) {
            argi.push_back(arg_split[j][i]);
          } else {
            argi.push_back(arg_split[j][0]);
          }
        }
        const_cast<Function&>(fcn)->call(argi, ret[i], false, false);
      }

      ret = swapIndices(ret);

      std::vector<MX> ret_cat;
      for (int i=0;i<fcn.nOut();++i) {
        ret_cat.push_back(horzcat(ret[i]));
      }
      return ret_cat;
    } else {
      // Get type of parallelization
      bool omp;
      if (parallelization.compare("openmp")==0) {
        omp = true;
      } else if (parallelization.compare("serial")==0) {
        omp = false;
      } else {
        casadi_error("Unsupported parallelization \"" << parallelization
                     << "\": Available options are expand|serial|openmp");
      }

      // Call the map
      std::vector<MX> v;
      if (omp) {
        return MX::createMultipleOutput(new OmpMap(fcn, arg));
      } else {
        return MX::createMultipleOutput(new Map(fcn, arg));
      }

    }
  }


} // namespace casadi
