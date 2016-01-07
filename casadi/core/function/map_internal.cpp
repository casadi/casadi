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
#include "../profiling.hpp"

using namespace std;

namespace casadi {

  MapBase* MapBase::create(const Function& f, int n, const Dict& opts) {
    // Check if there are reduced inputs or outputs
    bool reduce_inputs = opts.find("reduced_inputs")!=opts.end();
    bool reduce_outputs = opts.find("reduced_outputs")!=opts.end();

    // Read the type of parallelization
    Dict::const_iterator par_op = opts.find("parallelization");
    if (reduce_inputs || reduce_outputs) {
      // Vector indicating which inputs/outputs are to be repeated
      std::vector<bool> repeat_in(f.nIn(), true), repeat_out(f.nOut(), true);

      // Mark reduced inputs
      if (reduce_inputs) {
        vector<int> ri = opts.find("reduced_inputs")->second;
        for (vector<int>::const_iterator i=ri.begin(); i!=ri.end(); ++i) {
          repeat_in[*i]=false;
        }
      }

      // Mark reduced outputs
      if (reduce_inputs) {
        vector<int> ro = opts.find("reduced_outputs")->second;
        for (vector<int>::const_iterator i=ro.begin(); i!=ro.end(); ++i) {
          repeat_out[*i]=false;
        }
      }


      bool non_repeated_output = false;
      for (int i=0;i<repeat_out.size();++i) {
        non_repeated_output |= !repeat_out[i];
      }

      if (par_op==opts.end() || par_op->second == "serial") {
        return new MapSumSerial(f, n, repeat_in, repeat_out);
      } else {
        if (par_op->second == "openmp") {
          if (non_repeated_output) {
            casadi_warning("OpenMP not yet supported for reduced outputs. "
                           "Falling back to serial mode.");
          } else {
            #ifdef WITH_OPENMP
            return new MapSumOmp(f, n, repeat_in, repeat_out);
            #else // WITH_OPENMP
            casadi_warning("CasADi was not compiled with OpenMP. "
                           "Falling back to serial mode.");
            #endif // WITH_OPENMP
          }
        } else if (par_op->second == "opencl") {
          if (non_repeated_output) {
            casadi_warning("OpenCL not yet supported for reduced outputs. "
                           "Falling back to serial mode.");
          } else {
            return new MapSumOcl(f, n, repeat_in, repeat_out);
          }
        }
        return new MapSumSerial(f, n, repeat_in, repeat_out);
      }
    }

    if (par_op==opts.end() || par_op->second == "serial") {
      return new MapSerial(f, n);
    } else {
      if (par_op->second == "openmp") {
        #ifdef WITH_OPENMP
        return new MapOmp(f, n);
        #else // WITH_OPENMP
        casadi_warning("CasADi was not compiled with OpenMP. "
                       "Falling back to serial mode.");
        #endif // WITH_OPENMP
      } else if (par_op->second == "opencl") {
        return new MapOcl(f, n);
      }
      return new MapSerial(f, n);
    }
  }

  MapBase::MapBase(const Function& f, int n)
    : f_(f), n_in_(f.nIn()), n_out_(f.nOut()), n_(n) {

    addOption("parallelization", OT_STRING, "serial",
              "Computational strategy for parallelization", "serial|openmp|opencl");
    addOption("opencl_select", OT_INTEGER, 0,
      "List with indices into OpenCL-compatible devices, to select which one to use.");

    addOption("reduced_inputs", OT_INTEGERVECTOR, GenericType(),
              "Reduction for certain inputs");
    addOption("reduced_outputs", OT_INTEGERVECTOR, GenericType(),
              "Reduction for certain outputs");

    setOption("input_scheme", f.inputScheme());
    setOption("output_scheme", f.outputScheme());

    // Give a name
    setOption("name", "unnamed_map");

  }

  MapBase::~MapBase() {
  }

  void MapBase::init() {
    // Call the initialization method of the base class
    FunctionInternal::init();

    opencl_select_ = getOption("opencl_select");

  }

  void PureMap::init() {
    // Initialize the functions, get input and output sparsities
    // Input and output sparsities

    ibuf_.resize(f_.nIn());
    obuf_.resize(f_.nOut());

    for (int i=0;i<f_.nIn();++i) {
      // Sparsity of the original input
      Sparsity in_sp = f_.input(i).sparsity();

      // Allocate space for map input
      input(i) = DMatrix::zeros(repmat(in_sp, 1, n_));
    }

    for (int i=0;i<f_.nOut();++i) {
      // Sparsity of the original output
      Sparsity out_sp = f_.output(i).sparsity();

      // Allocate space for map output
      output(i) = DMatrix::zeros(repmat(out_sp, 1, n_));
    }

    MapBase::init();

    // Allocate sufficient memory for serial evaluation
    alloc_arg(f_.sz_arg());
    alloc_res(f_.sz_res());
    alloc_w(f_.sz_w());
    alloc_iw(f_.sz_iw());


    step_in_.resize(nIn(), 0);
    step_out_.resize(nOut(), 0);

    for (int i=0;i<nIn();++i) {
      step_in_[i] = f_.input(i).nnz();
    }

    for (int i=0;i<nOut();++i) {
      step_out_[i] = f_.output(i).nnz();
    }

  }

  template<typename T>
  void PureMap::evalGen(const T** arg, T** res, int* iw, T* w) {
    const T** arg1 = arg+n_in_;
    T** res1 = res+n_out_;

    for (int i=0; i<n_; ++i) {
      for (int j=0; j<n_in_;  ++j)
        arg1[j] = (arg[j]==0) ? 0: arg[j]+i*step_in_[j];
      for (int j=0; j<n_out_; ++j)
        res1[j] = (res[j]==0) ? 0: res[j]+i*step_out_[j];
      f_(arg1, res1, iw, w);
    }
  }

  void PureMap::evalSX(const SXElement** arg, SXElement** res, int* iw, SXElement* w) {
    evalGen(arg, res, iw, w);
  }

  void PureMap::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    evalGen(arg, res, iw, w);
  }

  void PureMap::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    bvec_t** arg1 = arg+n_in_;
    bvec_t** res1 = res+n_out_;

    for (int i=0; i<n_; ++i) {
      for (int j=0; j<n_in_;  ++j)
        arg1[j] = (arg[j]==0) ? 0: arg[j]+i*step_in_[j];
      for (int j=0; j<n_out_; ++j)
        res1[j] = (res[j]==0) ? 0: res[j]+i*step_out_[j];
      f_->spAdj(arg1, res1, iw, w);
    }
  }

  Function PureMap
  ::getDerForward(const std::string& name, int nfwd, Dict& opts) {

    // Differentiate mapped function
    Function df = f_.derForward(nfwd);

    // Propagate options (if not set already)
    if (opts.find("parallelization")==opts.end()) {
      opts["parallelization"] = parallelization();
    }
    if (opts.find("opencl_select")==opts.end()) {
      opts["opencl_select"] = opencl_select_;
    }

    // Construct and return
    return Map(name, df, n_, opts);
  }

  Function PureMap
  ::getDerReverse(const std::string& name, int nadj, Dict& opts) {
    // Differentiate mapped function
    Function df = f_.derReverse(nadj);

    // Propagate options (if not set already)
    if (opts.find("parallelization")==opts.end()) {
      opts["parallelization"] = parallelization();
    }
    if (opts.find("opencl_select")==opts.end()) {
      opts["opencl_select"] = opencl_select_;
    }

    // Construct and return
    return Map(name, df, n_, opts);
  }

  void PureMap::generateDeclarations(CodeGenerator& g) const {
    f_->addDependency(g);
  }

  void PureMap::generateBody(CodeGenerator& g) const {

    int num_in = f_.nIn(), num_out = f_.nOut();

    g.body << "  const real_t** arg1 = arg+" << f_.sz_arg() << ";" << endl;
    g.body << "  real_t** res1 = res+"  << f_.sz_res() <<  ";" << endl;

    g.body << "  int i;" << endl;
    g.body << "  for (i=0; i<" << n_ << "; ++i) {" << endl;

    for (int j=0; j<num_in; ++j) {
      g.body << "    arg1[" << j << "] = (arg[" << j << "]==0)? " <<
        "0: arg[" << j << "]+i*" << step_in_[j] << ";" << endl;
    }
    for (int j=0; j<num_out; ++j) {
      g.body << "    res1[" << j << "] = (res[" << j << "]==0)? " <<
        "0: res[" << j << "]+i*" << step_out_[j] << ";" << endl;
    }

    g.body << "    " << g.call(f_, "arg1", "res1", "iw", "w") << ";" << endl;
    g.body << "  }" << std::endl;
  }

  MapSum::MapSum(const Function& f, int n,
                       const std::vector<bool> &repeat_in,
                       const std::vector<bool> &repeat_out)
    : MapBase(f, n), repeat_in_(repeat_in), repeat_out_(repeat_out) {

    casadi_assert_message(repeat_in_.size()==f.nIn(),
                          "MapSum expected repeat_in of size " << f.nIn() <<
                          ", but got " << repeat_in_.size() << " instead.");

    casadi_assert_message(repeat_out_.size()==f.nOut(),
                          "MapSum expected repeat_out of size " << f.nOut() <<
                          ", but got " << repeat_out_.size() << " instead.");

  }

  MapSum::~MapSum() {

  }


  void MapSum::init() {

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
    MapBase::init();

    // Allocate some space to evaluate each function to.
    nnz_out_ = 0;
    for (int i=0;i<num_out;++i) {
      if (!repeat_out_[i]) nnz_out_+= step_out_[i];
    }

    alloc_w(f_.sz_w() + nnz_out_);
    alloc_iw(f_.sz_iw());
    alloc_arg(2*f_.sz_arg());
    alloc_res(2*f_.sz_res());

  }


  template<typename T, typename R>
  void MapSum::evalGen(const T** arg, T** res, int* iw, T* w,
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

  void MapSum::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    int num_in = f_.nIn(), num_out = f_.nOut();

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

      f_.spAdj(arg1, res1, iw, w);
    }

    // Reset all seeds
    for (int j=0; j<num_out; ++j) {
      if (res[j]!=0) {
        fill(res[j], res[j]+f_.output(j).nnz(), bvec_t(0));
      }
    }

  }


  void MapSum::evalSX(const SXElement** arg, SXElement** res,
                           int* iw, SXElement* w) {
    evalGen<SXElement>(arg, res, iw, w, &FunctionInternal::evalSX, std::plus<SXElement>());
  }

  static bvec_t Orring(bvec_t x, bvec_t y) { return x | y; }

  void MapSum::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    evalGen<bvec_t>(arg, res, iw, w, &FunctionInternal::spFwd, &Orring);
  }

  Function MapSum
  ::getDerForward(const std::string& name, int nfwd, Dict& opts) {

    // Differentiate mapped function
    Function df = f_.derForward(nfwd);

    // Propagate options (if not set already)
    if (opts.find("parallelization")==opts.end()) {
      opts["parallelization"] = parallelization();
    }
    if (opts.find("opencl_select")==opts.end()) {
      opts["opencl_select"] = opencl_select_;
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

  Function MapSum
  ::getDerReverse(const std::string& name, int nadj, Dict& opts) {
    // Differentiate mapped function
    Function df = f_.derReverse(nadj);

    // Propagate options (if not set already)
    if (opts.find("parallelization")==opts.end()) {
      opts["parallelization"] = parallelization();
    }
    if (opts.find("opencl_select")==opts.end()) {
      opts["opencl_select"] = opencl_select_;
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

  void MapSum::generateDeclarations(CodeGenerator& g) const {
    f_->addDependency(g);
  }

  void MapSum::generateBody(CodeGenerator& g) const {

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

  void MapSum::print(ostream &stream) const {
    stream << "Map(" << name(f_) << ", " << n_ << ")";
  }

  MapSerial::~MapSerial() {
  }

  void MapSerial::init() {
    // Call the initialization method of the base class
    PureMap::init();

  }

  void MapSerial::evalD(const double** arg, double** res, int* iw, double* w) {
    evalGen(arg, res, iw, w);
  }

  void MapSumSerial::init() {
    MapSum::init();
  }

  MapSumSerial::~MapSumSerial() {

  }

  void MapSumSerial::evalD(const double** arg, double** res,
                          int* iw, double* w) {
    double t0 = getRealTime();
    evalGen<double>(arg, res, iw, w, &FunctionInternal::eval, std::plus<double>());
  }

#ifdef WITH_OPENMP

  MapOmp::~MapOmp() {
  }

  void MapOmp::evalD(const double** arg, double** res, int* iw, double* w) {
    size_t sz_arg, sz_res, sz_iw, sz_w;
    f_.sz_work(sz_arg, sz_res, sz_iw, sz_w);
#pragma omp parallel for
    for (int i=0; i<n_; ++i) {
      int n_in = n_in_, n_out = n_out_;
      const double** arg_i = arg + n_in*n_ + sz_arg*i;
      copy(arg+i*n_in, arg+(i+1)*n_in, arg_i);
      double** res_i = res + n_out*n_ + sz_res*i;
      copy(res+i*n_out, res+(i+1)*n_out, res_i);
      int* iw_i = iw + i*sz_iw;
      double* w_i = w + i*sz_w;
      f_->evalD(arg_i, res_i, iw_i, w_i);
    }
  }

  void MapOmp::init() {
    // Call the initialization method of the base class
    PureMap::init();

    // Allocate sufficient memory for parallel evaluation
    alloc_arg(f_.sz_arg() * n_);
    alloc_res(f_.sz_res() * n_);
    alloc_w(f_.sz_w() * n_);
    alloc_iw(f_.sz_iw() * n_);
  }

  void MapSumOmp::evalD(const double** arg, double** res,
                          int* iw, double* w) {
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
        f_->eval(arg_i, res_i, iw_i, w_i);
      }
  }


  void MapSumOmp::init() {
    // Call the initialization method of the base class
    MapSum::init();

    // Allocate sufficient memory for parallel evaluation
    alloc_w(f_.sz_w()*n_ + nnz_out_);
    alloc_iw(f_.sz_iw()*n_);
    alloc_arg(f_.sz_arg()*(n_+1));
    alloc_res(f_.sz_res()*(n_+1));
  }

  MapSumOmp::~MapSumOmp() {

  }

#endif // WITH_OPENMP

  MapOcl::MapOcl(const Function& f, int n) : PureMap(f, n) {

  }

  MapOcl::~MapOcl() {

  }

  void MapOcl::generateDeclarations(CodeGenerator& g) const {
    f_->addDependency(g);
  }

  void MapOcl::generateBody(CodeGenerator& g) const {
    g.addAuxiliary(CodeGenerator::AUX_MINMAX);
    g.addAuxiliary(CodeGenerator::AUX_ASSERT);

    g.addInclude("CL/cl.h");

    // Get a unique index for thsi function
    int& ind = g.added_dependencies_[this];

    // Static global memory (will not be needed in CasADi 3.0)
    g.declarations << "static cl_kernel kernel" << ind << "_ = 0;" << endl;
    g.declarations << "static cl_command_queue commands" << ind << "_ = 0;" << endl;
    g.declarations << "static cl_context context" << ind << "_ = 0;" << endl;
    g.declarations << "static cl_program program" << ind << "_ = 0;" << endl;
    g.declarations << "static cl_mem d_in" << ind << "_ = 0;" << endl;
    g.declarations << "static cl_mem d_out" << ind << "_ = 0;" << endl;

    g.cleanup << "clReleaseMemObject(d_in"  << ind << "_);" << endl;
    g.cleanup << "clReleaseMemObject(d_out" << ind << "_);" << endl;
    g.cleanup << "clReleaseProgram(program" << ind << "_);" << endl;
    g.cleanup << "clReleaseKernel(kernel" << ind << "_);" << endl;
    g.cleanup << "clReleaseCommandQueue(commands" << ind << "_);" << endl;
    g.cleanup << "clReleaseContext(context" << ind << "_);" << endl;

    // Debug OpenCL errors
    g.declarations << "#define check_cl_error(a) ";
    g.declarations << "assert_action(a == CL_SUCCESS, printf(\"OpenCL exit code '%d'\\n\",a))";
    g.declarations << endl;

    // =========================
    // OpenCL setup code START
    // =========================

    g.setup << "  {" << endl;

    // Select the desirable OpenCL devices
    g.setup << "    cl_uint numPlatforms;" << endl;
    g.setup << "    int err = clGetPlatformIDs(0, NULL, &numPlatforms);" << endl;
    g.setup << "    check_cl_error(err);" << endl;

    g.setup << "    cl_platform_id Platform[numPlatforms];" << endl;
    g.setup << "    err = clGetPlatformIDs(numPlatforms, Platform, NULL);" << endl;
    g.setup << "    check_cl_error(err);" << endl;

    g.setup << "    cl_device_id mydevice;"  << endl;

    g.setup << "    int i,j,k=0;" << endl;
    g.setup << "    for (i = 0; i < numPlatforms; i++) {" << endl;
    g.setup << "      cl_uint n = 0;" << endl;
    g.setup << "      err = clGetDeviceIDs(Platform[i], CL_DEVICE_TYPE_ALL, 0, NULL, &n);" << endl;
    g.setup << "      check_cl_error(err);" << endl;
    g.setup << "      cl_device_id device_id[n];" << endl;
    g.setup << "      err = clGetDeviceIDs(Platform[i], CL_DEVICE_TYPE_ALL, n, device_id, NULL);";
    g.setup << endl;
    g.setup << "      check_cl_error(err);" << endl;
    g.setup << "      for (j=0;j<n;++j) {" << endl;
    g.setup << "        cl_char device_name[1024] = {0};" << endl;
    g.setup << "        err = clGetDeviceInfo(device_id[j], CL_DEVICE_NAME, sizeof(device_name),";
    g.setup << " &device_name, NULL);" << endl;
    g.setup << "        check_cl_error(err);" << endl;
    g.setup << "        printf(\"Detected device %d: %s\", k, device_name);" << endl;
    g.setup << "        if (k== " << opencl_select_ << ") {" << endl;
    g.setup << "          mydevice = device_id[j];" << std::endl;
    g.setup << "          printf(\" (selected)\");" << endl;
    g.setup << "        }" << endl;
    g.setup << "        printf(\"\\n\");" << endl;
    g.setup << "        k+=1;" << endl;
    g.setup << "      }" << endl;
    g.setup << "    }" << endl;


    // Create a compute context
    g.setup << "    context" << ind << "_ = clCreateContext(0, 1, &mydevice, NULL, NULL, &err);";
    g.setup << endl;
    g.setup << "    check_cl_error(err);" << endl;

    // Create a command queue
    g.setup << "    commands" << ind << "_ = clCreateCommandQueue(context" << ind << "_,";
    g.setup << " mydevice, 0, &err);" << endl;
    g.setup << "    check_cl_error(err);" << endl;

    // Obtain the kernel source
    g.setup << "    const char *KernelSource = " << g.multiline_string(kernelCode()+"\n") << ";";
    g.setup << endl;

    // Create the compute program from the source buffer
    g.setup << "    program" << ind << "_ = clCreateProgramWithSource(context" << ind << "_, 1,";
    g.setup << " (const char **) & KernelSource, NULL, &err);" << endl;


    // Build the program
    g.setup << "    err = clBuildProgram(program" << ind << "_, 0, NULL, NULL, NULL, NULL);";
    g.setup << endl;
    g.setup << "    if (err != CL_SUCCESS) {" << endl;
    g.setup << "      size_t len;" << endl;
    g.setup << "      char buffer[200048];" << endl;
    g.setup << "      clGetProgramBuildInfo(program" << ind << "_, mydevice, CL_PROGRAM_BUILD_LOG,";
    g.setup << " sizeof(buffer), buffer, &len);" << endl;
    g.setup << "      printf(\"%s\\n\", buffer);" << endl;
    g.setup << "      check_cl_error(err);" << endl;
    g.setup << "    }" << endl;


    // Create the compute kernel from the program
    g.setup << "    kernel" << ind << "_ = clCreateKernel(program" << ind;
    g.setup << "_, \"mykernel\", &err);" << endl;
    g.setup << "    check_cl_error(err);" << endl;

    // Create buffer objects
    g.setup << "    d_in" << ind << "_ = clCreateBuffer(context" << ind << "_, CL_MEM_READ_ONLY,";
    g.setup << " sizeof(float)*" <<  f_.nnzIn()*n_ << ", NULL, &err);" << endl;
    g.setup << "    check_cl_error(err);" << endl;
    g.setup << "    d_out" << ind << "_ = clCreateBuffer(context" << ind << "_, CL_MEM_WRITE_ONLY,";
    g.setup << " sizeof(float)*" << f_.nnzOut()*n_ << ", NULL, &err);" << endl;
    g.setup << "    check_cl_error(err);" << endl;

    // Prepare kernel arguments
    g.setup << "    err   = clSetKernelArg(kernel" << ind << "_, 0,";
    g.setup << " sizeof(cl_mem), &d_in" << ind << "_);" << endl;
    g.setup << "    err  |= clSetKernelArg(kernel" << ind << "_, 1,";
    g.setup << " sizeof(cl_mem), &d_out" << ind << "_);" << endl;


    g.setup << "  }" << endl;

    // =========================
    // OpenCL setup code END
    // =========================

    // =========================
    // OpenCL driver code START
    // =========================

    g.body << "  int i,j,k,kk;" << endl;
    g.body << "  size_t num = " << n_ << ";" << endl;

    // Obtain pointers to float host buffers
    g.body << "  float *h_in_ = (float*) (w+" << f_.sz_w() << ");" << endl;
    g.body << "  float *h_out_ = (float*) (w+" << f_.sz_w()  << "+sizeof(float)*(";
    g.body << f_.nnzIn()*n_ << ")/sizeof(double));" << endl;

    g.body << "  kk=0;" << endl;
    g.body << "  for (j=0;j<" << n_ << ";++j) {" << endl;
    for (int i=0;i<f_.nIn();++i) {
      g.body << "    if (arg[" << i << "]) {" << endl;
      g.body << "      for (k=0;k<" << f_.inputSparsity(i).nnz() << ";++k) {" << endl;
      g.body << "        h_in_[kk] = arg[" << i << "][k+j*" << f_.inputSparsity(i).nnz() << "];";
      g.body << endl;
      g.body << "        kk++;" << endl;
      g.body << "      }" << endl;
      g.body << "    } else {" << endl;
      g.body << "      k+= " << f_.inputSparsity(i).nnz() << ";" << endl;
      g.body << "    }" << endl;
    }
    g.body << "  }" << endl;

    // Copy data from host to buffer
    g.body << "  int err;" << std::endl;
    g.body << "  err = clEnqueueWriteBuffer(commands" << ind << "_, d_in" << ind << "_, CL_TRUE,";
    g.body << " 0, sizeof(float) * " << f_.nnzIn()*n_ << ",h_in_, 0, NULL, NULL);" << endl;
    g.body << "  check_cl_error(err);" << endl;

    // Run the kernel
    g.body << "  err = clEnqueueNDRangeKernel(commands" << ind << "_, kernel" << ind << "_, 1,";
    g.body << " NULL, &num, NULL, 0, NULL, NULL);" << endl;
    g.body << "  check_cl_error(err);" << endl;
    g.body << "  err = clFinish(commands" << ind << "_);" << endl;
    g.body << "  check_cl_error(err);" << endl;

    // Read back the output
    g.body << "  err = clEnqueueReadBuffer(commands" << ind << "_, d_out" << ind << "_, CL_TRUE,";
    g.body << " 0, " << sizeof(float) * (f_.nnzOut()*n_) << ", h_out_,";
    g.body << " 0, NULL, NULL ); " << endl;
    g.body << "  check_cl_error(err);" << endl;

    // Sum up the results
    g.body << "  kk=0;" << std::endl;
    g.body << "  for (j=0;j<" << n_ << ";++j) {" << std::endl;
    for (int i=0;i<f_.nOut();++i) {
      g.body << "    if (res[" << i << "]) {" << std::endl;
      g.body << "      for (k=0;k<" << f_.outputSparsity(i).nnz() << ";++k) {" << std::endl;
      g.body << "        res[" << i << "][k+j*" << f_.outputSparsity(i).nnz() << "] = h_out_[kk];";
      g.body << std::endl;
      g.body << "        kk++;" << std::endl;
      g.body << "      }" << std::endl;
      g.body << "    } else {" << std::endl;
      g.body << "      kk+= " << f_.outputSparsity(i).nnz() << ";" << std::endl;
      g.body << "    }" << std::endl;
    }
    g.body << "  }" << std::endl;

    // =========================
    // OpenCL driver code END
    // =========================

  }

  void MapOcl::evalD(const double** arg, double** res, int* iw, double* w) {
    casadi_error("OpenCL works only in codegeneration mode.");
  }

  std::string MapOcl::kernelCode() const {
    // Generate code for the kernel
    std::stringstream code;

    Dict opts;
    opts["opencl"] = true;
    opts["meta"] = false;

    CodeGenerator cg(opts);
    cg.add(f_, "F");

    //code << "#pragma OPENCL EXTENSION cl_khr_fp64: enable" << std::endl;
    code << "#define d float" << std::endl;
    code << "#define real_t float" << std::endl;
    code << "#define CASADI_PREFIX(ID) test_c_ ## ID" << std::endl;

    code << cg.generate() << std::endl;

    code << "__kernel void mykernel(" << std::endl;
    code << "   __global float* h_in," << std::endl;
    code << "   __global float* h_out)" << std::endl;
    code << "{                            " << std::endl;
    code << "   float h_in_local[" << f_.nnzIn() << "];" << std::endl;
    code << "   float h_out_local[" << f_.nnzOut() << "];" << std::endl;
    code << "   int i = get_global_id(0);" << std::endl;
    code << "   for (int k=0;k<"<< f_.nnzIn() << ";++k)" << std::endl;
    code << "     h_in_local[k] = h_in[k+i*" << f_.nnzIn() << "];" << std::endl;
    code << "   int iw[" << f_.sz_iw() << "];" << std::endl;
    code << "   float w[" << f_.sz_w() << "];" << std::endl;
    code << "   const d* arg[" << f_.sz_arg() << "];" << std::endl;
    code << "   d* res[" << f_.sz_res() << "];" << std::endl;

    int offset= 0;
    for (int i=0;i<f_.nIn();++i) {
      code << "   arg[" << i << "] = h_in_local+" << offset << ";" << std::endl;
      offset+= f_.inputSparsity(i).nnz();
    }
    offset= 0;
    for (int i=0;i<f_.nOut();++i) {
      code << "   res[" << i << "] = h_out_local+" << offset << ";" << std::endl;
      offset+= f_.outputSparsity(i).nnz();
    }

    code << "   F(arg,res,iw,w); " << std::endl;
    code << "   for (int k=0;k<"<< f_.nnzOut() << ";++k)" << std::endl;
    code << "     h_out[k+i*" << f_.nnzOut() << "] = h_out_local[k];" << std::endl;

    code << "}   " << std::endl;

    return code.str();
  }

  void MapOcl::init() {

    PureMap::init();

    alloc_w(f_.sz_w() +
      sizeof(float)*(f_.nnzIn()*n_ + f_.nnzOut()*n_)/sizeof(double));

  }

  MapSumOcl::MapSumOcl(const Function& f, int n,
    const std::vector<bool> &repeat_in, const std::vector<bool> &repeat_out) :
      MapSum(f, n, repeat_in, repeat_out) {

  }


  std::string MapSumOcl::kernelCode() const {
    // Generate code for the kernel
    std::stringstream code;

    Dict opts;
    opts["opencl"] = true;
    opts["meta"] = false;
    opts["main"] = true;

    CodeGenerator cg(opts);
    cg.add(f_, "F");

    //code << "#pragma OPENCL EXTENSION cl_khr_fp64: enable" << std::endl;
    code << "#define d float" << std::endl;
    code << "#define real_t float" << std::endl;
    code << "#define CASADI_PREFIX(ID) test_c_ ## ID" << std::endl;

    code << cg.generate() << std::endl;

    code << "__kernel void mykernel(" << std::endl;
    code << "   __global float* h_in," << std::endl;
    code << "   __global float* h_fixed," << std::endl;
    code << "   __global float* h_out)" << std::endl;
    code << "{                            " << std::endl;

    //code << "for (int i=0;i<" << f_.nnzOut()*n_ << ";++i) { h_out[i] = 1;}" << std::endl;
    code << "   float h_in_local[" << nnz_in_ << "];" << std::endl;
    code << "   float h_fixed_local[" << nnz_fixed_ << "];" << std::endl;
    code << "   float h_out_local[" << f_.nnzOut() << "];" << std::endl;
    code << "   int i = get_global_id(0);" << std::endl;
    code << "   for (int k=0;k<"<< nnz_in_ << ";++k)" << std::endl;
    code << "     h_in_local[k] = h_in[k+i*" << nnz_in_ << "];" << std::endl;
    code << "   for (int k=0;k<"<< nnz_fixed_ << ";++k)" << std::endl;
    code << "     h_fixed_local[k] = h_fixed[k];" << std::endl;
    code << "   int iw[" << f_.sz_iw()  << "];" << std::endl;
    code << "   float w[" << f_.sz_w() << "];" << std::endl;
    code << "   const d* arg[" << f_.sz_arg() << "];" << std::endl;
    code << "   d* res[" << f_.sz_res() << "];" << std::endl;

    int offset = 0;
    int offset_fixed = 0;
    for (int i=0;i<f_.nIn();++i) {
      if (repeat_in_[i]) {
        code << "   arg[" << i << "] = h_in_local+" << offset << ";" << std::endl;
        offset+= f_.inputSparsity(i).nnz();
      } else {
        code << "   arg[" << i << "] = h_fixed_local+" << offset_fixed << ";" << std::endl;
        offset_fixed += f_.inputSparsity(i).nnz();
      }
    }

    offset= 0;
    for (int i=0;i<f_.nOut();++i) {
      code << "   res[" << i << "] = h_out_local+" << offset  << ";" << std::endl;
      offset+= f_.outputSparsity(i).nnz();
    }

    code << "   F(arg, res, iw, w); " << std::endl;
    code << "   for (int k=0;k<"<< f_.nnzOut() << ";++k)" << std::endl;
    code << "     h_out[k+i*" << f_.nnzOut() << "] = h_out_local[k];" << std::endl;

    code << "}   " << std::endl;

    return code.str();
  }

  void MapSumOcl::evalD(const double** arg, double** res, int* iw, double* w) {
    casadi_error("OpenCL works only in codegeneration mode.");
  }

  void MapSumOcl::generateDeclarations(CodeGenerator& g) const {
    f_->addDependency(g);
  }

  void MapSumOcl::generateBody(CodeGenerator& g) const {
    g.addAuxiliary(CodeGenerator::AUX_MINMAX);
    g.addAuxiliary(CodeGenerator::AUX_ASSERT);

    g.addInclude("CL/cl.h");

    // Get a unique index for thsi function
    int& ind = g.added_dependencies_[this];

    // Static global memory (will not be needed in CasADi 3.0)
    g.declarations << "static cl_kernel kernel" << ind << "_ = 0;" << endl;
    g.declarations << "static cl_command_queue commands" << ind << "_ = 0;" << endl;
    g.declarations << "static cl_context context" << ind << "_ = 0;" << endl;
    g.declarations << "static cl_program program" << ind << "_ = 0;" << endl;
    g.declarations << "static cl_mem d_in" << ind << "_ = 0;" << endl;
    g.declarations << "static cl_mem d_fx" << ind << "_ = 0;" << endl;
    g.declarations << "static cl_mem d_out" << ind << "_ = 0;" << endl;

    g.cleanup << "clReleaseMemObject(d_in"  << ind << "_);" << endl;
    g.cleanup << "clReleaseMemObject(d_fx" << ind << "_);" << endl;
    g.cleanup << "clReleaseMemObject(d_out" << ind << "_);" << endl;
    g.cleanup << "clReleaseProgram(program" << ind << "_);" << endl;
    g.cleanup << "clReleaseKernel(kernel" << ind << "_);" << endl;
    g.cleanup << "clReleaseCommandQueue(commands" << ind << "_);" << endl;
    g.cleanup << "clReleaseContext(context" << ind << "_);" << endl;

    // Debug OpenCL errors
    g.declarations << "#define check_cl_error(a) ";
    g.declarations << "assert_action(a == CL_SUCCESS, printf(\"OpenCL exit code '%d'\\n\",a))";
    g.declarations << endl;

    // =========================
    // OpenCL setup code START
    // =========================

    g.setup << "  {" << endl;

    // Select the desirable OpenCL devices
    g.setup << "    cl_uint numPlatforms;" << endl;
    g.setup << "    int err = clGetPlatformIDs(0, NULL, &numPlatforms);" << endl;
    g.setup << "    check_cl_error(err);" << endl;

    g.setup << "    cl_platform_id Platform[numPlatforms];" << endl;
    g.setup << "    err = clGetPlatformIDs(numPlatforms, Platform, NULL);" << endl;
    g.setup << "    check_cl_error(err);" << endl;

    g.setup << "    cl_device_id mydevice;"  << endl;

    g.setup << "    int i,j,k=0;" << endl;
    g.setup << "    for (i = 0; i < numPlatforms; i++) {" << endl;
    g.setup << "      cl_uint n = 0;" << endl;
    g.setup << "      err = clGetDeviceIDs(Platform[i], CL_DEVICE_TYPE_ALL, 0, NULL, &n);" << endl;
    g.setup << "      check_cl_error(err);" << endl;
    g.setup << "      cl_device_id device_id[n];" << endl;
    g.setup << "      err = clGetDeviceIDs(Platform[i], CL_DEVICE_TYPE_ALL, n, device_id, NULL);";
    g.setup << endl;
    g.setup << "      check_cl_error(err);" << endl;
    g.setup << "      for (j=0;j<n;++j) {" << endl;
    g.setup << "        cl_char device_name[1024] = {0};" << endl;
    g.setup << "        err = clGetDeviceInfo(device_id[j], CL_DEVICE_NAME, sizeof(device_name),";
    g.setup << " &device_name, NULL);" << endl;
    g.setup << "        check_cl_error(err);" << endl;
    g.setup << "        printf(\"Detected device %d: %s\", k, device_name);" << endl;
    g.setup << "        if (k== " << opencl_select_ << ") {" << endl;
    g.setup << "          mydevice = device_id[j];" << std::endl;
    g.setup << "          printf(\" (selected)\");" << endl;
    g.setup << "        }" << endl;
    g.setup << "        printf(\"\\n\");" << endl;
    g.setup << "        k+=1;" << endl;
    g.setup << "      }" << endl;
    g.setup << "    }" << endl;


    // Create a compute context
    g.setup << "    context" << ind << "_ = clCreateContext(0, 1, &mydevice, NULL, NULL, &err);";
    g.setup << endl;
    g.setup << "    check_cl_error(err);" << endl;

    // Create a command queue
    g.setup << "    commands" << ind << "_ = clCreateCommandQueue(context" << ind << "_,";
    g.setup << " mydevice, 0, &err);" << endl;
    g.setup << "    check_cl_error(err);" << endl;

    // Obtain the kernel source
    g.setup << "    const char *KernelSource = " << g.multiline_string(kernelCode()+"\n") << ";";
    g.setup << endl;

    // Create the compute program from the source buffer
    g.setup << "    program" << ind << "_ = clCreateProgramWithSource(context" << ind << "_, 1,";
    g.setup << " (const char **) & KernelSource, NULL, &err);" << endl;


    // Build the program
    g.setup << "    err = clBuildProgram(program" << ind << "_, 0, NULL, NULL, NULL, NULL);";
    g.setup << endl;
    g.setup << "    if (err != CL_SUCCESS) {" << endl;
    g.setup << "      size_t len;" << endl;
    g.setup << "      char buffer[200048];" << endl;
    g.setup << "      clGetProgramBuildInfo(program" << ind << "_, mydevice, CL_PROGRAM_BUILD_LOG,";
    g.setup << " sizeof(buffer), buffer, &len);" << endl;
    g.setup << "      printf(\"%s\\n\", buffer);" << endl;
    g.setup << "      check_cl_error(err);" << endl;
    g.setup << "    }" << endl;


    // Create the compute kernel from the program
    g.setup << "    kernel" << ind << "_ = clCreateKernel(program" << ind;
    g.setup << "_, \"mykernel\", &err);" << endl;
    g.setup << "    check_cl_error(err);" << endl;

    // Create buffer objects
    g.setup << "    d_in" << ind << "_ = clCreateBuffer(context" << ind << "_, CL_MEM_READ_ONLY,";
    g.setup << " sizeof(float)*" <<  nnz_in_*n_ << ", NULL, &err);" << endl;
    g.setup << "    check_cl_error(err);" << endl;
    g.setup << "    d_fx" << ind << "_ = clCreateBuffer(context" << ind << "_, CL_MEM_READ_ONLY,";
    g.setup << " sizeof(float)*" <<  nnz_fixed_ << ", NULL, &err);" << endl;
    g.setup << "    check_cl_error(err);" << endl;

    g.setup << "    d_out" << ind << "_ = clCreateBuffer(context" << ind << "_, CL_MEM_WRITE_ONLY,";
    g.setup << " sizeof(float)*" << f_.nnzOut()*n_ << ", NULL, &err);" << endl;
    g.setup << "    check_cl_error(err);" << endl;

    // Prepare kernel arguments
    g.setup << "    err   = clSetKernelArg(kernel" << ind << "_, 0,";
    g.setup << " sizeof(cl_mem), &d_in" << ind << "_);" << endl;
    g.setup << "    err   = clSetKernelArg(kernel" << ind << "_, 1,";
    g.setup << " sizeof(cl_mem), &d_fx" << ind << "_);" << endl;
    g.setup << "    err  |= clSetKernelArg(kernel" << ind << "_, 2,";
    g.setup << " sizeof(cl_mem), &d_out" << ind << "_);" << endl;


    g.setup << "  }" << endl;

    // =========================
    // OpenCL setup code END
    // =========================

    // =========================
    // OpenCL driver code START
    // =========================

    g.body << "  int i,j,k,kk;" << endl;
    g.body << "  size_t num = " << n_ << ";" << endl;

    // Obtain pointers to float host buffers
    g.body << "  float *h_in_ = (float*) (w+" << f_.sz_w() << ");" << endl;
    g.body << "  float *h_fx_ = (float*) (w+" << f_.sz_w()  << "+sizeof(float)*(";
    g.body << nnz_in_*n_ << ")/sizeof(double));" << endl;
    g.body << "  float *h_out_ = (float*) (w+" << f_.sz_w()  << "+sizeof(float)*(";
    g.body << nnz_in_*n_ +  nnz_fixed_  << ")/sizeof(double));" << endl;

    g.body << "  kk=0;" << endl;
    g.body << "  for (j=0;j<" << n_ << ";++j) {" << endl;
    for (int i=0;i<f_.nIn();++i) {
      if (repeat_in_[i]) {
        g.body << "    if (arg[" << i << "]) {" << endl;
        g.body << "      for (k=0;k<" << f_.inputSparsity(i).nnz() << ";++k) {" << endl;
        g.body << "        h_in_[kk] = arg[" << i << "][k+j*" << f_.inputSparsity(i).nnz() << "];";
        g.body << endl;
        g.body << "        kk++;" << endl;
        g.body << "      }" << endl;
        g.body << "    } else {" << endl;
        g.body << "      k+= " << f_.inputSparsity(i).nnz() << ";" << endl;
        g.body << "    }" << endl;
      }
    }
    g.body << "  }" << endl;

    for (int i=0;i<f_.nIn();++i) {
      if (!repeat_in_[i]) {
        g.body << "  if (arg[" << i << "]) {" << endl;
        g.body << "    for (k=0;k<" << f_.inputSparsity(i).nnz() << ";++k) {" << endl;
        g.body << "      h_fx_[kk] = arg[" << i << "][k+j*" << f_.inputSparsity(i).nnz() << "];";
        g.body << endl;
        g.body << "      kk++;" << endl;
        g.body << "    }" << endl;
        g.body << "  } else {" << endl;
        g.body << "    k+= " << f_.inputSparsity(i).nnz() << ";" << endl;
        g.body << "  }" << endl;
      }
    }

    // Copy data from host to buffer
    g.body << "  int err;" << std::endl;
    g.body << "  err = clEnqueueWriteBuffer(commands" << ind << "_, d_in" << ind << "_, CL_TRUE,";
    g.body << " 0, sizeof(float) * " << nnz_in_*n_ << ",h_in_, 0, NULL, NULL);" << endl;
    g.body << "  check_cl_error(err);" << endl;
    g.body << "  err = clEnqueueWriteBuffer(commands" << ind << "_, d_fx" << ind << "_, CL_TRUE,";
    g.body << " 0, sizeof(float) * " << nnz_fixed_ << ",h_fx_, 0, NULL, NULL);" << endl;
    g.body << "  check_cl_error(err);" << endl;

    // Run the kernel
    g.body << "  err = clEnqueueNDRangeKernel(commands" << ind << "_, kernel" << ind << "_, 1,";
    g.body << " NULL, &num, NULL, 0, NULL, NULL);" << endl;
    g.body << "  check_cl_error(err);" << endl;
    g.body << "  err = clFinish(commands" << ind << "_);" << endl;
    g.body << "  check_cl_error(err);" << endl;

    // Read back the output
    g.body << "  err = clEnqueueReadBuffer(commands" << ind << "_, d_out" << ind << "_, CL_TRUE,";
    g.body << " 0, sizeof(float) *" <<  (f_.nnzOut()*n_) << ", h_out_,";
    g.body << " 0, NULL, NULL ); " << endl;
    g.body << "  check_cl_error(err);" << endl;

    // Sum up the results
    g.body << "  kk=0;" << std::endl;
    g.body << "  for (j=0;j<" << n_ << ";++j) {" << std::endl;
    for (int i=0;i<f_.nOut();++i) {
      g.body << "    if (res[" << i << "]) {" << std::endl;
      g.body << "      for (k=0;k<" << f_.outputSparsity(i).nnz() << ";++k) {" << std::endl;
      g.body << "        res[" << i << "][k+j*" << f_.outputSparsity(i).nnz() << "] = h_out_[kk];";
      g.body << std::endl;
      g.body << "        kk++;" << std::endl;
      g.body << "      }" << std::endl;
      g.body << "    } else {" << std::endl;
      g.body << "      kk+= " << f_.outputSparsity(i).nnz() << ";" << std::endl;
      g.body << "    }" << std::endl;
    }
    g.body << "  }" << std::endl;

    // =========================
    // OpenCL driver code END
    // =========================

  }

  void MapSumOcl::init() {
    MapSum::init();

    nnz_in_ = 0;nnz_fixed_ = 0;
    for (int i=0;i<n_in_;++i) {
      if (repeat_in_[i]) {
        nnz_in_  += f_.input(i).nnz();
      } else {
        nnz_fixed_ += f_.input(i).nnz();
      }
    }

    if (nnz_fixed_ == 0) nnz_fixed_ = 1;

    alloc_w(f_.sz_w() +
      sizeof(float)*(nnz_in_*n_ + nnz_fixed_ + f_.nnzOut()*n_)/sizeof(double));

  }

  MapSumOcl::~MapSumOcl() {

  }

} // namespace casadi
