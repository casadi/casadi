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


#include "map.hpp"

using namespace std;

namespace casadi {

  Function MapBase::create(const std::string& name,
                          const std::string& parallelization, Function& f, int n,
                          const std::vector<int>& reduce_in, const std::vector<int>& reduce_out,
                          const Dict& opts) {
    // Sanity checks on arguments
    casadi_assert(n>0);
    casadi_assert(inBounds(reduce_in, f.n_in()) && isUnique(reduce_in));
    casadi_assert(inBounds(reduce_out, f.n_out()) && isUnique(reduce_out));

    if (parallelization=="unroll") {
      // Obtain lookup tables
      std::vector<int> reduce_in_compl = complement(reduce_in, f.n_in());
      std::vector<int> reduce_out_compl = complement(reduce_out, f.n_out());

      // Construct symbolic inputs
      std::vector< std::vector<MX> > args(f.n_in());
      std::vector< MX > args_local(f.n_in());
      for (int k : reduce_in) args_local[k] = MX::sym("x", f.sparsity_in(k));
      for (int k : reduce_in_compl) {
        for (int i=0;i<n;++i) args[k].reserve(n);
        for (int i=0;i<n;++i) args[k].push_back(MX::sym("x", f.sparsity_in(k)));
      }

      // Collect results of symbolic calls to f
      std::vector< std::vector<MX> > ret(f.n_out());
      for (int k=0;k<f.n_out();++k) ret[k].reserve(n);
      for (int i=0;i<n;++i) {
        for (int k : reduce_in_compl) args_local[k] = args[k][i];
        std::vector<MX> ret_local = f(args_local);
        for (int k=0;k<f.n_out();++k) ret[k].push_back(ret_local[k]);
      }

      // Produce an MX function as a result
      std::vector< MX > f_in(f.n_in());
      for (int k : reduce_in) f_in[k] = args_local[k];
      for (int k : reduce_in_compl) f_in[k] = horzcat(args[k]);
      std::vector< MX > f_out(f.n_out());
      for (int k : reduce_out) f_out[k] = repsum(horzcat(ret[k]), 1, n);
      for (int k : reduce_out_compl) f_out[k] = horzcat(ret[k]);

      return Function(name, f_in, f_out, opts);
    } else {
      Function ret;
      ret.assignNode(MapBase::create(name, parallelization, f, n, reduce_in, reduce_out));
      ret->construct(opts);
      return ret;
    }
  }

  MapBase* MapBase::create(const std::string& name,
                          const std::string& parallelization, const Function& f, int n,
                          const std::vector<int>& reduce_in, const std::vector<int>& reduce_out) {

    if (reduce_in.size()>0 || reduce_out.size()>0) {
      // Vector indicating which inputs/outputs are to be repeated
      std::vector<bool> repeat_in(f.n_in(), true), repeat_out(f.n_out(), true);
      for (int i : reduce_in) repeat_in[i]= false;
      for (int i : reduce_out) repeat_out[i]= false;
      return new MapSum(name, f, n, repeat_in, repeat_out);
    } else {
      if (parallelization == "serial") {
        return new Map(name, f, n);
      } else if (parallelization== "openmp") {
        return new MapOmp(name, f, n);
      } else {
        casadi_error("Unknown parallelization: " + parallelization);
      }
    }
  }

  MapBase::MapBase(const std::string& name, const Function& f, int n)
    : FunctionInternal(name), f_(f), n_in_(f.n_in()), n_out_(f.n_out()), n_(n) {
  }

  Options MapBase::options_
  = {{&FunctionInternal::options_},
     {{"n_threads",
       {OT_INT,
        "Control the number of threads when executing in parallel. "
        "The default setting (0) will pass the decision on to the parallelization library. "
        "For openmp, this means that OMP_NUM_THREADS env. variable is observed."}}
     }
  };

  void MapBase::propagate_options(Dict& opts) {
    if (opts.find("n_threads")==opts.end()) opts["n_threads"] = n_threads_;
  }

  MapBase::~MapBase() {
  }

  void Map::init(const Dict& opts) {
    // Call the initialization method of the base class
    MapBase::init(opts);

    // Allocate sufficient memory for serial evaluation
    alloc_arg(f_.sz_arg());
    alloc_res(f_.sz_res());
    alloc_w(f_.sz_w());
    alloc_iw(f_.sz_iw());
  }

  template<typename T>
  void Map::evalGen(const T** arg, T** res, int* iw, T* w) const {
    int n_in = n_in_, n_out = n_out_;
    const T** arg1 = arg+this->n_in();
    T** res1 = res+this->n_out();
    for (int i=0; i<n_; ++i) {
      for (int j=0; j<n_in; ++j) {
        arg1[j] = arg[j] ? arg[j]+i*f_.nnz_in(j): 0;
      }
      for (int j=0; j<n_out; ++j) {
        res1[j]= res[j] ? res[j]+i*f_.nnz_out(j): 0;
      }
      f_(arg1, res1, iw, w, 0);
    }
  }

  void Map::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    evalGen(arg, res, iw, w);
  }

  void Map::sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    evalGen(arg, res, iw, w);
  }

  void Map::sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    int n_in = n_in_, n_out = n_out_;
    bvec_t** arg1 = arg+this->n_in();
    bvec_t** res1 = res+this->n_out();
    for (int i=0; i<n_; ++i) {
      for (int j=0; j<n_in; ++j) {
        arg1[j] = arg[j] ? arg[j]+i*f_.nnz_in(j): 0;
      }
      for (int j=0; j<n_out; ++j) {
        res1[j]= res[j] ? res[j]+i*f_.nnz_out(j): 0;
      }
      f_->sp_rev(arg1, res1, iw, w, 0);
    }
  }

  void Map::generateDeclarations(CodeGenerator& g) const {
    f_->addDependency(g);
  }

  void Map::generateBody(CodeGenerator& g) const {

    g.body << "  const real_t** arg1 = arg+" << n_in() << ";"<< endl;
    g.body << "  real_t** res1 = res+" << n_out() << ";" << endl;

    g.body << "  int i;" << endl;
    g.body << "  for (i=0; i<" << n_ << "; ++i) {" << endl;
    for (int j=0; j<n_in_; ++j) {
      g.body << "    arg1[" << j << "] = arg[" << j << "]? " <<
        "arg[" << j << "]+i*" << f_.nnz_in(j) << " : 0;" << endl;
    }
    for (int j=0; j<n_out_; ++j) {
      g.body << "    res1[" << j << "] = res[" << j << "]? " <<
        "res[" << j << "]+i*" << f_.nnz_out(j) << " : 0;" << endl;
    }
    g.body << "    if (" << g(f_, "arg1", "res1", "iw", "w") << ") return 1;" << endl;
    g.body << "  }" << std::endl;
  }

  Function Map
  ::get_forward_old(const std::string& name, int nfwd, Dict& opts) {
    // Differentiate mapped function
    Function df = f_.forward(nfwd);

    // Propagate options
    propagate_options(opts);

    // Construct and return
    return df.map(name, parallelization(), n_, opts);
  }

  Function Map
  ::get_reverse_old(const std::string& name, int nadj, Dict& opts) {
    // Differentiate mapped function
    Function df = f_.reverse(nadj);

    // Propagate options
    propagate_options(opts);

    // Construct and return
    return df.map(name, parallelization(), n_, opts);
  }

  MapSum::MapSum(const std::string& name, const Function& f, int n,
                       const std::vector<bool> &repeat_in,
                       const std::vector<bool> &repeat_out)
    : MapBase(name, f, n), repeat_in_(repeat_in), repeat_out_(repeat_out) {

    casadi_assert_message(repeat_in_.size()==f.n_in(),
                          "MapSum expected repeat_in of size " << f.n_in() <<
                          ", but got " << repeat_in_.size() << " instead.");

    casadi_assert_message(repeat_out_.size()==f.n_out(),
                          "MapSum expected repeat_out of size " << f.n_out() <<
                          ", but got " << repeat_out_.size() << " instead.");

  }

  MapSum::~MapSum() {

  }

  void MapBase::init(const Dict& opts) {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    // Read the 'n_threads' option and do a sanity check
    if (opts.find("n_threads")!=opts.end()) {
      n_threads_ = opts.find("n_threads")->second;
    } else {
      n_threads_ = 0;
    }
    casadi_assert_message(n_threads_>=0, "'n_threads' option must be a positive integer.");

  }

  void MapSum::init(const Dict& opts) {

    int num_in = f_.n_in(), num_out = f_.n_out();

    // Initialize the functions, get input and output sparsities
    // Input and output sparsities
    step_in_.resize(num_in, 0);
    step_out_.resize(num_out, 0);
    for (int i=0; i<num_in;++i) {
      if (repeat_in_[i]) step_in_[i] = f_.nnz_in(i);
    }
    for (int i=0; i<num_out; ++i) step_out_[i] = f_.nnz_out(i);

    // Call the initialization method of the base class
    MapBase::init(opts);

    // Allocate some space to evaluate each function to.
    nnz_out_ = 0;
    for (int i=0; i<num_out; ++i) {
      if (!repeat_out_[i]) nnz_out_+= step_out_[i];
    }

    alloc_w(f_.sz_w() + nnz_out_);
    alloc_iw(f_.sz_iw());
    alloc_arg(2*f_.sz_arg());
    alloc_res(2*f_.sz_res());
  }

  template<typename T, typename R>
  void MapSum::evalGen(const T** arg, T** res, int* iw, T* w, R reduction) const {
    int num_in = f_.n_in(), num_out = f_.n_out();

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
        arg1[j] = arg[j] ? arg[j]+i*step_in_[j] : 0;
      }

      // Set the function outputs
      for (int j=0; j<num_out; ++j) {
        if (repeat_out_[j]) {
          // Make the function outputs end up in our outputs
          res1[j] = res[j] ? res[j]+i*step_out_[j]: 0;
        } else {
          // Make the function outputs end up in temp_res
          res1[j] = res[j] ? temp_res : 0;
          temp_res+= step_out_[j];
        }
      }

      // Evaluate the function
      f_(arg1, res1, iw, w, 0);

      // Sum results from temporary storage to accumulator
      for (int k=0;k<num_out;++k) {
        if (res1[k] && sum[k] && !repeat_out_[k])
          std::transform(res1[k], res1[k]+step_out_[k], sum[k], sum[k], reduction);
      }
    }
  }

  void MapSum::eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    evalGen<double>(arg, res, iw, w, std::plus<double>());
  }

  void MapSum::eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
    evalGen<SXElem>(arg, res, iw, w, std::plus<SXElem>());
  }

  static bvec_t Orring(bvec_t x, bvec_t y) { return x | y; }

  void MapSum::sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    evalGen<bvec_t>(arg, res, iw, w, &Orring);
  }

  void MapSum::sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    int num_in = f_.n_in(), num_out = f_.n_out();

    bvec_t** arg1 = arg+f_.sz_arg();
    bvec_t** res1 = res+f_.sz_res();

    for (int i=0; i<n_; ++i) {

      bvec_t* temp_res = w+f_.sz_w();

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
          if (res[j]!=0) {
            copy(res[j], res[j]+step_out_[j], temp_res);
          }
          temp_res+= step_out_[j];
        }
      }

      f_->sp_rev(arg1, res1, iw, w, 0);
    }

    // Reset all seeds
    for (int j=0; j<num_out; ++j) {
      if (res[j]!=0) {
        fill(res[j], res[j]+f_.nnz_out(j), bvec_t(0));
      }
    }

  }

  Function MapSum
  ::get_forward(const std::string& name, int nfwd, Dict& opts) {
    return get_forward_new(name, nfwd, opts);
  }

  Function MapSum
  ::get_forward_old(const std::string& name, int nfwd, Dict& opts) {

    // Differentiate mapped function
    Function df = f_.forward(nfwd);

    // Propagate options
    propagate_options(opts);

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

    std::vector<int> reduce_in;
    for (int i=0;i<repeat_in.size();++i) {
      if (!repeat_in[i]) reduce_in.push_back(i);
    }
    std::vector<int> reduce_out;
    for (int i=0;i<repeat_out.size();++i) {
      if (!repeat_out[i]) reduce_out.push_back(i);
    }

    // Construct and return
    return df.map(name, parallelization(), n_, reduce_in, reduce_out, opts);
  }

  Function MapSum
  ::get_reverse(const std::string& name, int nadj, Dict& opts) {
    return get_reverse_new(name, nadj, opts);
  }

  Function MapSum
  ::get_reverse_old(const std::string& name, int nadj, Dict& opts) {
    // Differentiate mapped function
    Function df = f_.reverse(nadj);

    // Propagate options
    propagate_options(opts);

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

    std::vector<int> reduce_in;
    for (int i=0;i<repeat_in.size();++i) {
      if (!repeat_in[i]) reduce_in.push_back(i);
    }
    std::vector<int> reduce_out;
    for (int i=0;i<repeat_out.size();++i) {
      if (!repeat_out[i]) reduce_out.push_back(i);
    }

    // Construct and return
    return df.map(name, parallelization(), n_, reduce_in, reduce_out, opts);
  }

  void MapSum::generateDeclarations(CodeGenerator& g) const {
    f_->addDependency(g);
  }

  void MapSum::generateBody(CodeGenerator& g) const {

    int num_in = f_.n_in(), num_out = f_.n_out();

    g.body << "  const real_t** arg1 = arg+" << f_.sz_arg() << ";" << endl
           << "  real_t** sum = res;" << endl;

    // Clear the accumulators
    for (int k=0;k<num_out;++k) {
      g.body << "  if (sum[" << k << "]!=0) " <<
        g.fill(STRING("sum[" << k << "]"), step_out_[k], "0") << endl;
    }

    g.body << "  real_t** res1 = res+"  << f_.sz_res() <<  ";" << endl;

    g.body << "  int i;" << endl;
    g.body << "  for (i=0; i<" << n_ << "; ++i) {" << endl;

    g.body << "    real_t* temp_res = w+"  << f_.sz_w() <<  ";" << endl
           << "    if (temp_res!=0) " << g.fill("temp_res", nnz_out_, "0") << endl;

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

    g.body << "    " << g(f_, "arg1", "res1", "iw", "w") << ";" << endl;

    g.addAuxiliary(CodeGenerator::AUX_AXPY);
    // Sum results
    for (int k=0; k<num_out; ++k) {
      if (!repeat_out_[k]) {
        g.body << "    if (res1[" << k << "] && sum[" << k << "])" << endl
               << "       axpy(" << step_out_[k] << ",1," <<
          "res1["<< k << "],sum[" << k << "]);" << endl;
      }
    }
    g.body << "  }" << std::endl;
  }

  void MapSum::print(ostream &stream) const {
    stream << "Map(" << f_.name() << ", " << n_ << ")";
  }

  void Map::eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    evalGen(arg, res, iw, w);
  }

  MapOmp::~MapOmp() {
  }

  void MapOmp::eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
#ifdef WITH_OPENMP
    return Map::eval(mem, arg, res, iw, w);
#else // WITH_OPENMP
    size_t sz_arg, sz_res, sz_iw, sz_w;
    f_.sz_work(sz_arg, sz_res, sz_iw, sz_w);

    if (n_threads_==0) {
#pragma omp parallel for
      for (int i=0; i<n_; ++i) {
        const double** arg_i = arg + n_in_ + sz_arg*i;
        for (int j=0; j<n_in_; ++j) {
          arg_i[j] = arg[j]+i*f_.nnz_in(j);
        }
        double** res_i = res + n_out_ + sz_res*i;
        for (int j=0; j<n_out_; ++j) {
          res_i[j] = res[j]? res[j]+i*f_.nnz_out(j) : 0;
        }
        int* iw_i = iw + i*sz_iw;
        double* w_i = w + i*sz_w;
        f_->eval(0, arg_i, res_i, iw_i, w_i);
      }
    } else {
#pragma omp parallel for num_threads(n_threads_)
      for (int i=0; i<n_; ++i) {
        const double** arg_i = arg + n_in_ + sz_arg*i;
        for (int j=0; j<n_in_; ++j) {
          arg_i[j] = arg[j]+i*f_.nnz_in(j);
        }
        double** res_i = res + n_out_ + sz_res*i;
        for (int j=0; j<n_out_; ++j) {
          res_i[j] = res[j]? res[j]+i*f_.nnz_out(j) : 0;
        }
        int* iw_i = iw + i*sz_iw;
        double* w_i = w + i*sz_w;
        f_->eval(0, arg_i, res_i, iw_i, w_i);
      }
    }
  #endif  // WITH_OPENMP
  }

  void MapOmp::generateDeclarations(CodeGenerator& g) const {
    f_->addDependency(g);
  }

  void MapOmp::generateBody(CodeGenerator& g) const {
    size_t sz_arg, sz_res, sz_iw, sz_w;
    f_.sz_work(sz_arg, sz_res, sz_iw, sz_w);

    g.body << "  int i;" << endl;
    if (n_threads_==0) {
      g.body << "#pragma omp parallel for" << endl;
    } else {
      g.body << "#pragma omp parallel for num_threads(" << n_threads_ << ")" << endl;
    }
    g.body << "  for (i=0; i<" << n_ << "; ++i) {" << endl;
    g.body << "    const double** arg_i = arg + " << n_in_ << "+" << sz_arg << "*i;" << endl;
    for (int j=0; j<n_in_; ++j) {
      g.body << "    arg_i[" << j << "] = arg[" << j << "]+i*" << f_.nnz_in(j) << ";" << endl;
    }
    g.body << "    double** res_i = res + " <<  n_out_ << "+" <<  sz_res << "*i;" << endl;
    for (int j=0; j<n_out_; ++j) {
      g.body << "    res_i[" << j << "] = res[" << j << "] ?" <<
                "res[" << j << "]+i*" << f_.nnz_out(j) << ": 0;" << endl;
    }
    g.body << "    int* iw_i = iw + i*" << sz_iw << ";" << endl;
    g.body << "    double* w_i = w + i*" << sz_w << ";" << endl;
    g.body << "    " << g(f_, "arg_i", "res_i", "iw_i", "w_i") << ";" << endl;
    g.body << "  }" << std::endl;
  }

  void MapOmp::init(const Dict& opts) {
    // Call the initialization method of the base class
    Map::init(opts);

    // Allocate sufficient memory for parallel evaluation
    alloc_arg(f_.sz_arg() * n_);
    alloc_res(f_.sz_res() * n_);
    alloc_w(f_.sz_w() * n_);
    alloc_iw(f_.sz_iw() * n_);
  }

} // namespace casadi
