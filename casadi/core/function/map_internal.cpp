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


#include "map_internal.hpp"
#include "mx_function.hpp"

using namespace std;

namespace casadi {

  MapInternal::MapInternal(const Function& f, int n,
    const std::vector<bool> &repeat_in,
    const std::vector<bool> &repeat_out)
    : f_(f), n_(n), repeat_in_(repeat_in), repeat_out_(repeat_out) {

    addOption("parallelization", OT_STRING, "serial",
              "Computational strategy for parallelization", "serial|openmp");

    casadi_assert_message(repeat_in_.size()==f.nIn(),
      "MapInternal expected repeat_in of size " << f.nIn() <<
      ", but got " << repeat_in_.size() << " instead.");

    casadi_assert_message(repeat_out_.size()==f.nOut(),
      "MapInternal expected repeat_out of size " << f.nOut() <<
      ", but got " << repeat_out_.size() << " instead.");

    // Give a name
    setOption("name", "unnamed_map");
  }

  MapInternal::~MapInternal() {
  }

  void MapInternal::init() {

    std::string parallelization = getOption("parallelization").toString();

    if (parallelization.compare("serial")==0) {
      parallelization_ = PARALLELIZATION_SERIAL;
    } else if (parallelization.compare("openmp")==0) {
      parallelization_ = PARALLELIZATION_OMP;
    } else {
      casadi_error("cannot happen.");
    }

    // OpenMP not yet supported for non-repeated outputs
    bool non_repeated_output = false;
    for (int i=0;i<repeat_out_.size();++i) {
      non_repeated_output |= !repeat_out_[i];
    }

    if (parallelization_ == PARALLELIZATION_OMP && non_repeated_output) {
      casadi_warning("OpenMP not yet supported for non-repeated outputs. " <<
                     "Falling back to serial mode.");
      parallelization_ = PARALLELIZATION_SERIAL;
    }
    #ifndef WITH_OPENMP
    if (parallelization_ == PARALLELIZATION_OMP) {
      casadi_warning("CasADi was not compiled with OpenMP." <<
                     "Falling back to serial mode.");
      parallelization_ = PARALLELIZATION_SERIAL;
    }
    #endif // WITH_OPENMP

    int num_in = f_.nIn(), num_out = f_.nOut();

    // Initialize the functions, get input and output sparsities
    // Input and output sparsities

    ibuf_.resize(num_in);
    obuf_.resize(num_out);

    for (int i=0;i<num_in;++i) {
      // Sparsity of the original input
      Sparsity in_sp = f_.input(i).sparsity();

      // Repeat if requested
      if (repeat_in_[i]) in_sp = repmat(in_sp, 1, n_);

      // Allocate space for input
      input(i) = DMatrix::zeros(in_sp);
    }

    for (int i=0;i<num_out;++i) {
      // Sparsity of the original output
      Sparsity out_sp = f_.output(i).sparsity();

      // Repeat if requested
      if (repeat_out_[i]) out_sp = repmat(out_sp, 1, n_);

      // Allocate space for output
      output(i) = DMatrix::zeros(out_sp);
    }

    step_in_.resize(num_in, 0);
    step_out_.resize(num_out, 0);

    for (int i=0;i<num_in;++i) {
        if (repeat_in_[i])
          step_in_[i] = f_.input(i).nnz();
    }

    for (int i=0;i<num_out;++i) {
      step_out_[i] = f_.output(i).nnz();
    }

    // Call the initialization method of the base class
    FunctionInternal::init();

    // Allocate some space to evaluate each function to.
    nnz_out_ = 0;
    for (int i=0;i<num_out;++i) {
      if (!repeat_out_[i]) nnz_out_+= step_out_[i];
    }

    if (parallelization_ == PARALLELIZATION_SERIAL) {
      alloc_w(f_.sz_w() + nnz_out_);
      alloc_iw(f_.sz_iw());
      alloc_arg(2*f_.sz_arg());
      alloc_res(2*f_.sz_res());
    } else if (parallelization_ == PARALLELIZATION_OMP) {
      alloc_w(f_.sz_w()*n_);
      alloc_iw(f_.sz_iw()*n_);
      alloc_arg(f_.sz_arg()*(n_+1));
      alloc_res(f_.sz_res()*(n_+1));
    }
  }

  template<typename T, typename R>
  void MapInternal::evalGen(const T** arg, T** res, int* iw, T* w,
    void (FunctionInternal::*eval)(const T** arg, T** res, int* iw, T* w),
    R reduction) {
    int num_in = f_.nIn(), num_out = f_.nOut();

    const T** arg1 = arg+f_.sz_arg();

    // Clear the accumulators
    T** sum = res;
    for (int k=0;k<num_out;++k) {
      if (sum[k]!=0) std::fill(sum[k], sum[k]+step_out_[k], 0);
    }

    T** res1 = res+f_.sz_res();

    for (int i=0; i<n_; ++i) {

      T* temp_res = w+f_.sz_w();
      // Clear the temp_res storage space
      if (temp_res!=0) std::fill(temp_res, temp_res+nnz_out_, 0);

      // Set the function inputs
      for (int j=0; j<num_in; ++j) {
        arg1[j] = (arg[j]==0) ? 0: arg[j]+i*step_in_[j];
      }

      // Set the function outputs
      for (int j=0; j<num_out; ++j) {
        if (repeat_out_[j]) {
          // Make the function outputs end up in our outputs
          res1[j] = (res[j]==0)? 0: res[j]+i*step_out_[j];
        } else {
          // Make the function outputs end up in temp_res
          res1[j] = (res[j]==0)? 0: temp_res;
          temp_res+= step_out_[j];
        }
      }

      // Evaluate the function
      FunctionInternal *f = f_.operator->();
      (f->*eval)(arg1, res1, iw, w);

      // Sum results from temporary storage to accumulator
      for (int k=0;k<num_out;++k) {
        if (res1[k] && sum[k] && !repeat_out_[k])
          std::transform(res1[k], res1[k]+step_out_[k], sum[k], sum[k], reduction);
      }
    }
  }

  void MapInternal::evalD(const double** arg, double** res,
                                int* iw, double* w) {
    if (parallelization_ == PARALLELIZATION_SERIAL) {
      evalGen<double>(arg, res, iw, w, &FunctionInternal::evalD, std::plus<double>());
    } else {
      int n_in_ = f_.nIn(), n_out_ = f_.nOut();
      #ifndef WITH_OPENMP
          casadi_error("the \"impossible\" happened: " <<
                       "should have fallen back to serial in init.");
      #endif // WITH_OPENMP
      size_t sz_arg, sz_res, sz_iw, sz_w;
      f_.sz_work(sz_arg, sz_res, sz_iw, sz_w);
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
        f_->evalD(arg_i, res_i, iw_i, w_i);
      }
    }

  }

  void MapInternal::evalSX(const SXElement** arg, SXElement** res,
                                int* iw, SXElement* w) {
    evalGen<SXElement>(arg, res, iw, w, &FunctionInternal::evalSX, std::plus<SXElement>());
  }

  static bvec_t Orring(bvec_t x, bvec_t y) { return x | y; }

  void MapInternal::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    evalGen<bvec_t>(arg, res, iw, w, &FunctionInternal::spFwd, &Orring);
  }

/**void Map::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    evalGen<SXElement>(arg, res, iw, w);
}*/

  /**
  void MapInternal::evalMX(const std::vector<MX>& arg, std::vector<MX>& res) {
    int num_in = f_.nIn(), num_out = f_.nOut();

      // split arguments
      std::vector< std::vector< MX > > arg_split;
      for (int i=0;i<arg.size();++i) {
        arg_split.push_back(horzsplit(arg[i], f_.input(i).size2()));
      }

      std::vector< std::vector<MX> > ret(n_);

      // call the original function
      for (int i=0; i<n_; ++i) {
        std::vector< MX > argi;
        for (int j=0;j<arg.size();++j) {
          if (repeat_in_[j]) {
            argi.push_back(arg_split[j][i]);
          } else {
            argi.push_back(arg_split[j][0]);
          }
        }
        const_cast<Function&>(f_)->call(argi, ret[i], false, false);
      }

      ret = transpose(ret);

      res.resize(num_out);

      std::vector<MX> ret_cat;
      for (int i=0;i<num_out;++i) {
        if (repeat_out_[i]) {
          res[i] = horzcat(ret[i]);
        } else {
          MX sum = 0;
          for (int j=0;j<ret[i].size();++j) {
            sum += ret[i][j];
          }
          res[i] = sum;
        }
      }

  }
*/

  Function MapInternal
  ::getDerForward(const std::string& name, int nfwd, Dict& opts) {

    // Differentiate mapped function
    Function df = f_.derForward(nfwd);

    // Propagate options (if not set already)
    if (opts.find("parallelization")==opts.end()) {
      switch (parallelization_) {
      case PARALLELIZATION_SERIAL: opts["parallelization"] = "serial"; break;
      case PARALLELIZATION_OMP: opts["parallelization"] = "openmp"; break;
      }
    }

    std::vector<bool> repeat_in;
    repeat_in.insert(repeat_in.end(), repeat_in_.begin(), repeat_in_.end());
    repeat_in.insert(repeat_in.end(), repeat_out_.begin(), repeat_out_.end());
    for (int i=0;i<nfwd;++i) {
      repeat_in.insert(repeat_in.end(), repeat_in_.begin(), repeat_in_.end());
    }

    std::vector<bool> repeat_out;
    for (int i=0;i<nfwd;++i) {
      repeat_out.insert(repeat_out.end(), repeat_out_.begin(), repeat_out_.end());
    }

    // Construct and return
    return Map(name, df, n_, repeat_in, repeat_out, opts);
  }

  Function MapInternal
  ::getDerReverse(const std::string& name, int nadj, Dict& opts) {
    // Differentiate mapped function
    Function df = f_.derReverse(nadj);

    // Propagate options (if not set already)
    if (opts.find("parallelization")==opts.end()) {
      switch (parallelization_) {
      case PARALLELIZATION_SERIAL: opts["parallelization"] = "serial"; break;
      case PARALLELIZATION_OMP: opts["parallelization"] = "openmp"; break;
      }
    }

    std::vector<bool> repeat_in;
    repeat_in.insert(repeat_in.end(), repeat_in_.begin(), repeat_in_.end());
    repeat_in.insert(repeat_in.end(), repeat_out_.begin(), repeat_out_.end());
    for (int i=0; i<nadj; ++i) {
      repeat_in.insert(repeat_in.end(), repeat_out_.begin(), repeat_out_.end());
    }

    std::vector<bool> repeat_out;
    for (int i=0; i<nadj; ++i) {
      repeat_out.insert(repeat_out.end(), repeat_in_.begin(), repeat_in_.end());
    }

    // Construct and return
    return Map(name, df, n_, repeat_in, repeat_out, opts);
  }

  void MapInternal::generateDeclarations(CodeGenerator& g) const {
    f_->addDependency(g);
  }

  void MapInternal::generateBody(CodeGenerator& g) const {

    int num_in = f_.nIn(), num_out = f_.nOut();

    g.body << "  const real_t** arg1 = arg+" << f_.sz_arg() << ";" << endl
           << "  real_t** sum = res;" << endl;

    // Clear the accumulators
    for (int k=0;k<num_out;++k) {
      g.body << "  if (sum[" << k << "]!=0) " <<
        g.fill_n(STRING("sum[" << k << "]"), step_out_[k], "0") << endl;
    }

    g.body << "  real_t** res1 = res+"  << f_.sz_res() <<  ";" << endl;

    g.body << "  int i;" << endl;
    g.body << "  for (i=0; i<" << n_ << "; ++i) {" << endl;

    g.body << "    real_t* temp_res = w+"  << f_.sz_w() <<  ";" << endl
           << "    if (temp_res!=0)" << g.fill_n("temp_res", nnz_out_, "0") << endl;

    for (int j=0; j<num_in; ++j) {
      g.body << "    arg1[" << j << "] = (arg[" << j << "]==0)? " <<
                "0: arg[" << j << "]+i*" << step_in_[j] << ";" << endl;
    }
    for (int j=0; j<num_out; ++j) {
      if (repeat_out_[j]) {
        g.body << "    res1[" << j << "] = (res[" << j << "]==0)? " <<
         "0: res[" << j << "]+i*" << step_out_[j] << ";" << endl;
      } else {
        g.body << "    res1[" << j << "] = (res[" << j << "]==0)? 0: temp_res;" << endl
               << "    temp_res+= " << step_out_[j] << ";" << endl;
      }
    }

    g.body << "    " << g.call(f_, "arg1", "res1", "iw", "w") << ";" << endl;

    g.addAuxiliary(CodeGenerator::AUX_AXPY);
    // Sum results
    for (int k=0; k<num_out; ++k) {
      if (!repeat_out_[k]) {
        g.body << "    if (res1[" << k << "] && sum[" << k << "])" << endl
               << "       axpy(" << step_out_[k] << ",1," <<
                            "res1["<< k << "],1,sum[" << k << "],1);" << endl;
      }
    }
    g.body << "  }" << std::endl;
  }

  inline string name(const Function& f) {
    if (f.isNull()) {
      return "NULL";
    } else {
      return f.getOption("name");
    }
  }

  void MapInternal::print(ostream &stream) const {
    stream << "Map(" << name(f_) << ", " << n_ << ")";
  }

} // namespace casadi
