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


#include "function_internal.hpp"
#include "../std_vector_tools.hpp"
#include "sx_function.hpp"
#include "mx_function.hpp"
#include "map.hpp"
#include "mapaccum.hpp"
#include "switch.hpp"
#include "kernel_sum.hpp"
#include "nlpsol.hpp"
#include "conic.hpp"
#include "jit.hpp"
#include "../casadi_file.hpp"

#include <typeinfo>
#include <fstream>
#include <cctype>

using namespace std;

namespace casadi {

  Function::Function() {
  }

  Function::~Function() {
  }

  bool Function::proceed_to(std::istream& file, const std::string& str) {
    // Make sure that the file is ready for reading
    if (!file.good()) return false;
    // Have we already wrapped around once?
    //bool wrapped_around = false;
    // Read line-by-line
    string tmp;
    while (true) {
      // Read a word
      streampos cur_pos = file.tellg();
      file >> tmp;
      if (!file.good()) return false;

      // Check if match
      if (str==tmp) return true;

      // If comment, continue to the end of the line
      if (tmp.at(0)=='#') {
        file.ignore(numeric_limits<streamsize>::max(), '\n');
        continue;
      }

      // If mismatching name, rewind and break
      file.seekg(cur_pos);
      return false;
    }
  }

  Function::Function(const std::string& fname) {
    // Parse the file
    ParsedFile file(fname);

    // Create the corresponding class
    string classname = file.to_string("CLASS");

    if (classname=="Jit") {
      *this = jit(file);
    } else {
      casadi_error("Unknown Function type: " + classname);
    }
  }

  Function::Function(const string& name,
                     const std::vector<SX>& arg, const std::vector<SX>& res,
                     const Dict& opts) {
    construct(name, arg, res, opts);
  }

  Function::Function(const string& name,
                     const std::vector<SX>& arg, const std::vector<SX>& res,
                     const std::vector<string>& argn, const std::vector<string>& resn,
                     const Dict& opts) {
    construct(name, arg, res, argn, resn, opts);
  }

  Function::Function(const string& name,
                     const std::vector<MX>& arg, const std::vector<MX>& res,
                     const Dict& opts) {
    construct(name, arg, res, opts);
  }

  Function::Function(const string& name,
                     const std::vector<MX>& arg, const std::vector<MX>& res,
                     const std::vector<string>& argn, const std::vector<string>& resn,
                     const Dict& opts) {
    construct(name, arg, res, argn, resn, opts);
  }

  Function::Function(const string& name, SXIList arg, const SXVector& res, const Dict& opts) {
    construct(name, SXVector(arg), res, opts);
  }

  Function::Function(const string& name, const SXVector& arg, SXIList res, const Dict& opts) {
    construct(name, arg, SXVector(res), opts);
  }

  Function::Function(const string& name, SXIList arg, SXIList res, const Dict& opts) {
    construct(name, SXVector(arg), SXVector(res), opts);
  }

  Function::Function(const string& name, SXIList arg, const SXVector& res,
                     const StringVector& argn, const StringVector& resn, const Dict& opts) {
    construct(name, SXVector(arg), res, argn, resn, opts);
  }

  Function::Function(const string& name, const SXVector& arg, SXIList res,
                     const StringVector& argn, const StringVector& resn, const Dict& opts) {
    construct(name, arg, SXVector(res), argn, resn, opts);
  }

  Function::Function(const string& name, SXIList arg, SXIList res,
                     const StringVector& argn, const StringVector& resn, const Dict& opts) {
    construct(name, SXVector(arg), SXVector(res), argn, resn, opts);
  }

  Function::Function(const string& name, MXIList arg, const MXVector& res, const Dict& opts) {
    construct(name, MXVector(arg), res, opts);
  }

  Function::Function(const string& name, const MXVector& arg, MXIList res, const Dict& opts) {
    construct(name, arg, MXVector(res), opts);
  }

  Function::Function(const string& name, MXIList arg, MXIList res, const Dict& opts) {
    construct(name, MXVector(arg), MXVector(res), opts);
  }

  Function::Function(const string& name, MXIList arg, const MXVector& res,
                     const StringVector& argn, const StringVector& resn, const Dict& opts) {
    construct(name, MXVector(arg), res, argn, resn, opts);
  }

  Function::Function(const string& name, const MXVector& arg, MXIList res,
                     const StringVector& argn, const StringVector& resn, const Dict& opts) {
    construct(name, arg, MXVector(res), argn, resn, opts);
  }

  Function::Function(const string& name, MXIList arg, MXIList res,
                     const StringVector& argn, const StringVector& resn, const Dict& opts) {
    construct(name, MXVector(arg), MXVector(res), argn, resn, opts);
  }

  Function::Function(const string& name, const std::map<string, SX>& dict,
                     const vector<string>& argn, const vector<string>& resn,
                     const Dict& opts) {
    construct(name, dict, argn, resn, opts);
  }

  Function::Function(const string& name, const std::map<string, MX>& dict,
                     const vector<string>& argn, const vector<string>& resn,
                     const Dict& opts) {
    construct(name, dict, argn, resn, opts);
  }

  template<typename M>
  void Function::construct(const string& name, const std::map<string, M>& dict,
                           const vector<string>& argn,
                           const vector<string>& resn,
                           const Dict& opts) {
    vector<M> arg(argn.size()), res(resn.size());
    for (auto&& i : dict) {
      vector<string>::const_iterator it;
      if ((it=find(argn.begin(), argn.end(), i.first))!=argn.end()) {
        // Input expression
        arg[it-argn.begin()] = i.second;
      } else if ((it=find(resn.begin(), resn.end(), i.first))!=resn.end()) {
        // Output expression
        res[it-resn.begin()] = i.second;
      } else {
        // Neither
        casadi_error("Unknown dictionary entry: '" + i.first + "'");
      }
    }
    construct(name, arg, res, argn, resn, opts);
  }

  void Function::construct(const string& name,
                           const vector<SX>& arg, const vector<SX>& res,
                           const Dict& opts) {
    assignNode(new SXFunction(name, arg, res));
    (*this)->construct(opts);
  }

  void Function::construct(const string& name,
                           const vector<MX>& arg, const vector<MX>& res,
                           const Dict& opts) {
    assignNode(new MXFunction(name, arg, res));
    (*this)->construct(opts);
  }

  template<typename M>
  void Function::construct(const string& name,
                           const vector<M>& arg, const vector<M>& res,
                           const vector<string>& argn, const vector<string>& resn,
                           const Dict& opts) {
    Dict opts2 = opts;
    opts2["input_scheme"] = argn;
    opts2["output_scheme"] = resn;
    construct(name, arg, res, opts2);
  }

  Function Function::expand() const {
    return expand(name());
  }

  Function Function::expand(const string& name, const Dict& opts) const {
    vector<SX> arg = sx_in();
    vector<SX> res = Function(*this)(arg);
    vector<string> name_in = this->name_in();
    vector<string> name_out = this->name_out();
    Dict opts2(opts);
    if (!name_in.empty() && !opts.count("input_scheme")) opts2["input_scheme"]=name_in;
    if (!name_out.empty() && !opts.count("output_scheme")) opts2["output_scheme"]=name_out;
    return Function(name, arg, res, opts2);
  }

  Function Function::create(FunctionInternal* node) {
    Function ret;
    ret.assignNode(node);
    return ret;
  }

  FunctionInternal* Function::operator->() const {
    return get();
  }

  FunctionInternal* Function::get() const {
    return static_cast<FunctionInternal*>(SharedObject::get());
  }

  void Function::call(const vector<DM> &arg, vector<DM> &res,
                      bool always_inline, bool never_inline) {
    (*this)->call(arg, res, always_inline, never_inline);
  }

  void Function::call(const vector<SX> &arg, vector<SX>& res,
                      bool always_inline, bool never_inline) {
    (*this)->call(arg, res, always_inline, never_inline);
  }

  void Function::call(const vector<MX> &arg, vector<MX>& res,
                      bool always_inline, bool never_inline) {
    (*this)->call(arg, res, always_inline, never_inline);
  }

  vector<const double*> Function::buf_in(Function::VecArg arg) const {
    casadi_assert(arg.size()==n_in());
    auto arg_it=arg.begin();
    vector<const double*> buf_arg(sz_arg());
    for (unsigned int i=0; i<arg.size(); ++i) {
      casadi_assert(arg_it->size()==nnz_in(i));
      buf_arg[i] = get_ptr(*arg_it++);
    }
    return buf_arg;
  }

  vector<double*> Function::buf_out(Function::VecRes res) const {
    res.resize(n_out());
    auto res_it=res.begin();
    vector<double*> buf_res(sz_res());
    for (unsigned int i=0; i<res.size(); ++i) {
      res_it->resize(nnz_out(i));
      buf_res[i] = get_ptr(*res_it++);
    }
    return buf_res;
  }

  vector<double*> Function::buf_out(Function::VPrRes res) const {
    casadi_assert(res.size()==n_out());
    auto res_it=res.begin();
    vector<double*> buf_res(sz_res());
    for (unsigned int i=0; i<res.size(); ++i) {
      casadi_assert(*res_it!=0);
      (*res_it)->resize(nnz_out(i));
      buf_res[i] = get_ptr(**res_it++);
    }
    return buf_res;
  }

  vector<const double*> Function::buf_in(Function::MapArg arg) const {
    // Return value (RVO)
    vector<const double*> ret(sz_arg(), 0);

    // Read inputs
    for (auto i=arg.begin(); i!=arg.end(); ++i) {
      int ind = index_in(i->first);
      casadi_assert(i->second.size()==nnz_in(ind));
      ret[ind] = get_ptr(i->second);
    }

    return ret;
  }

  vector<double*> Function::buf_out(Function::MapRes res) const {
    // Return value (RVO)
    vector<double*> ret(sz_res(), 0);

    // Read outputs
    for (auto i=res.begin(); i!=res.end(); ++i) {
      int ind = index_out(i->first);
      i->second.resize(nnz_out(ind));
      ret[ind] = get_ptr(i->second);
    }

    return ret;
  }

  vector<double*> Function::buf_out(Function::MPrRes res) const {
    // Return value (RVO)
    vector<double*> ret(sz_res(), 0);

    // Read outputs
    for (auto i=res.begin(); i!=res.end(); ++i) {
      int ind = index_out(i->first);
      casadi_assert(i->second!=0);
      i->second->resize(nnz_out(ind));
      ret[ind] = get_ptr(*i->second);
    }

    return ret;
  }

  template<typename D>
  void Function::_call(vector<const D*> arg, vector<D*> res) {
    // Input buffer
    casadi_assert(arg.size()>=n_in());
    arg.resize(sz_arg());

    // Output buffer
    casadi_assert(res.size()>=n_out());
    res.resize(sz_res());

    // Work vectors
    vector<int> iw(sz_iw());
    vector<D> w(sz_w());

    // Evaluate memoryless
    (*this)(get_ptr(arg), get_ptr(res), get_ptr(iw), get_ptr(w), 0);
  }


  void Function::operator()(vector<const double*> arg, vector<double*> res) {
    return _call(arg, res);
  }

  void Function::operator()(vector<const bvec_t*> arg, vector<bvec_t*> res) {
    return _call(arg, res);
  }

  void Function::operator()(vector<const SXElem*> arg, vector<SXElem*> res) {
    return _call(arg, res);
  }

  void Function::rev(std::vector<bvec_t*> arg, std::vector<bvec_t*> res) {
    // Input buffer
    casadi_assert(arg.size()>=n_in());
    arg.resize(sz_arg());

    // Output buffer
    casadi_assert(res.size()>=n_out());
    res.resize(sz_res());

    // Work vectors
    vector<int> iw(sz_iw());
    vector<bvec_t> w(sz_w());

    // Evaluate memoryless
    rev(get_ptr(arg), get_ptr(res), get_ptr(iw), get_ptr(w), 0);
  }

  Function Function::mapaccum(const string& name, int n, const Dict& opts) {
    std::vector<int> accum_in = std::vector<int>(1, 0);
    std::vector<int> accum_out = std::vector<int>(1, 0);
    return mapaccum(name, n, accum_in, accum_out, opts);
  }

  Function Function::mapaccum(const string& name, int n,
                              const vector<int>& accum_in,
                              const vector<int>& accum_out,
                              const Dict& opts) {
    // Inherit the scheme
    Dict opts2 = opts;
    if (opts.find("input_scheme")==opts.end()) opts2["input_scheme"] = name_in();
    if (opts.find("output_scheme")==opts.end()) opts2["output_scheme"] = name_out();

    return Mapaccum::create(name, *this, n, accum_in, accum_out, opts2);
  }

  Function Function::mapaccum(const string& name, int n,
                              const vector<std::string>& accum_in,
                              const vector<std::string>& accum_out,
                              const Dict& opts) {
    std::vector<int> accum_in_num;
    for (int i=0;i<accum_in.size();++i) accum_in_num.push_back(index_in(accum_in[i]));
    std::vector<int> accum_out_num;
    for (int i=0;i<accum_out.size();++i) accum_out_num.push_back(index_out(accum_out[i]));
    return mapaccum(name, n, accum_in_num, accum_out_num, opts);
  }

  Function Function::map(const string& name, const std::string& parallelization, int n,
      const std::vector<int>& reduce_in, const std::vector<int>& reduce_out,
      const Dict& opts) {
    return MapBase::create(name, parallelization, *this, n, reduce_in, reduce_out, opts);
  }

  Function Function::map(const string& name, const std::string& parallelization, int n,
      const std::vector<std::string>& reduce_in, const std::vector<std::string>& reduce_out,
      const Dict& opts) {
    std::vector<int> reduce_in_num;
    for (int i=0;i<reduce_in.size();++i) reduce_in_num.push_back(index_in(reduce_in[i]));
    std::vector<int> reduce_out_num;
    for (int i=0;i<reduce_out.size();++i) reduce_out_num.push_back(index_out(reduce_out[i]));
    return map(name, parallelization, n, reduce_in_num, reduce_out_num, opts);
  }

  Function Function::map(const string& name, const std::string& parallelization, int n,
      const Dict& opts) {
    std::vector<int> reduce_in;
    std::vector<int> reduce_out;
    return MapBase::create(name, parallelization, *this, n, reduce_in, reduce_out, opts);
  }

  Function Function::slice(const std::vector<int>& order_in, const std::vector<int>& order_out,
      const Dict& opts) {
    // Get symbolic inputs
    std::vector<MX> in = mx_in();

    // Get symbolic outputs
    std::vector<MX> out = (*this)(in);

    // Shuffle the inputs
    std::vector<MX> new_in;
    for (int k : order_in) new_in.push_back(in.at(k));

    // Shuffle the outputs
    std::vector<MX> new_out;
    for (int k : order_out) new_out.push_back(out.at(k));

    // Return a wrapping function
    return Function(std::string("slice_") + name(), new_in, new_out, opts);
  }

  vector<MX> Function::map(const vector< MX > &x,
                           const string& parallelization) {
    return (*this)->map_mx(x, parallelization);
  }


  std::map<std::string, MX> Function::map(const std::map<std::string, MX> &arg,
                      const std::string& parallelization) {

    // Convert inputs map to vector
    std::vector<MX> args(n_in());
    for (auto& e : arg) args[index_in(e.first)] = e.second;

    // Delegate the actual map call
    std::vector<MX> res = map(args, parallelization);

    // Result vector to map
    std::map<std::string, MX> ret;
    for (int i=0;i<res.size();++i) ret[name_out(i)] = res[i];

    return ret;
  }

  vector<MX> Function::mapsum(const vector< MX > &x,
                              const string& parallelization) {
    return (*this)->mapsum_mx(x, parallelization);
  }

  Function Function::conditional(const string& name, const vector<Function>& f,
                                 const Function& f_def, const Dict& opts) {
    Function ret;
    ret.assignNode(new Switch(name, f, f_def));
    ret->construct(opts);
    return ret;
  }

  Function Function::if_else(const string& name, const Function& f_true,
                             const Function& f_false, const Dict& opts) {
    Function ret;
    ret.assignNode(new Switch(name, vector<Function>(1, f_false), f_true));
    ret->construct(opts);
    return ret;
  }

  Function Function::kernel_sum(const string& name,
                                const pair<int, int> & size,
                                double r, int n,
                                const Dict& opts) const {
    Function ret;
    ret.assignNode(new KernelSum(name, *this, size, r, n));
    ret->construct(opts);
    return ret;
  }

  int Function::n_in() const {
    return (*this)->n_in();
  }

  int Function::n_out() const {
    return (*this)->n_out();
  }

  int Function::size1_in(int ind) const {
    return (*this)->size1_in(ind);
  }

  int Function::size2_in(int ind) const {
    return (*this)->size2_in(ind);
  }

  int Function::size1_out(int ind) const {
    return (*this)->size1_out(ind);
  }

  int Function::size2_out(int ind) const {
    return (*this)->size2_out(ind);
  }

  pair<int, int> Function::size_in(int ind) const {
    return (*this)->size_in(ind);
  }

  pair<int, int> Function::size_out(int ind) const {
    return (*this)->size_out(ind);
  }

  int Function::nnz_in() const {
    return (*this)->nnz_in();
  }

  int Function::nnz_out() const {
    return (*this)->nnz_out();
  }

  int Function::numel_in() const {
    return (*this)->numel_in();
  }

  int Function::numel_out() const {
    return (*this)->numel_out();
  }

  int Function::nnz_in(int ind) const {
    return (*this)->nnz_in(ind);
  }

  int Function::nnz_out(int ind) const {
    return (*this)->nnz_out(ind);
  }

  int Function::numel_in(int ind) const {
    return (*this)->numel_in(ind);
  }

  int Function::numel_out(int ind) const {
    return (*this)->numel_out(ind);
  }

  Function Function::jacobian(int iind, int oind, bool compact, bool symmetric) {
    return (*this)->jacobian(iind, oind, compact, symmetric);
  }

  void Function::setJacobian(const Function& jac, int iind, int oind, bool compact) {
    (*this)->setJacobian(jac, iind, oind, compact);
  }

  Function Function::gradient(int iind, int oind) {
    return (*this)->gradient(iind, oind);
  }

  Function Function::tangent(int iind, int oind) {
    return (*this)->tangent(iind, oind);
  }

  Function Function::hessian(int iind, int oind) {
    return (*this)->hessian(iind, oind);
  }

  Function Function::fullJacobian() {
    return (*this)->fullJacobian();
  }

  void Function::setFullJacobian(const Function& jac) {
    (*this)->full_jacobian_ = jac;
  }

  bool Function::test_cast(const SharedObjectNode* ptr) {
    return dynamic_cast<const FunctionInternal*>(ptr)!=0;
  }

  Dict Function::stats(int mem) const {
    return (*this)->_get_stats(mem);
  }

  const Sparsity Function::sparsity_jac(int iind, int oind, bool compact, bool symmetric) {
    return (*this)->sparsity_jac(iind, oind, compact, symmetric);
  }

  void Function::set_jac_sparsity(const Sparsity& sp, int iind, int oind, bool compact) {
    (*this)->set_jac_sparsity(sp, iind, oind, compact);
  }

  vector<string> Function::name_in() const {
    return (*this)->ischeme_;
  }

  vector<string> Function::name_out() const {
    return (*this)->oscheme_;
  }

  int Function::index_in(const string &name) const {
    return (*this)->index_in(name);
  }

  int Function::index_out(const string &name) const {
    return (*this)->index_out(name);
  }

  string Function::name_in(int ind) const {
    return (*this)->name_in(ind);
  }

  string Function::name_out(int ind) const {
    return (*this)->name_out(ind);
  }

  const Sparsity& Function::sparsity_in(int ind) const {
    return (*this)->sparsity_in(ind);
  }

  const Sparsity& Function::sparsity_in(const string &iname) const {
    return (*this)->sparsity_in(iname);
  }

  const Sparsity& Function::sparsity_out(int ind) const {
    return (*this)->sparsity_out(ind);
  }

  const Sparsity& Function::sparsity_out(const string &iname) const {
    return (*this)->sparsity_out(iname);
  }

  void Function::sz_work(size_t& sz_arg, size_t& sz_res, size_t& sz_iw, size_t& sz_w) const {
    (*this)->sz_work(sz_arg, sz_res, sz_iw, sz_w);
  }

  size_t Function::sz_arg() const { return (*this)->sz_arg();}

  size_t Function::sz_res() const { return (*this)->sz_res();}

  size_t Function::sz_iw() const { return (*this)->sz_iw();}

  size_t Function::sz_w() const { return (*this)->sz_w();}

  void Function::operator()(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const {
    (*const_cast<Function*>(this))->sp_fwd(arg, res, iw, w, mem);
  }

  void Function::rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    (*this)->sp_rev(arg, res, iw, w, mem);
  }

  void Function::set_work(const double**& arg, double**& res, int*& iw, double*& w,
                          int mem) const {
    (*this)->_set_work(arg, res, iw, w, mem);
  }

  void Function::set_temp(const double** arg, double** res, int* iw, double* w,
                          int mem) const {
    (*this)->_set_temp(arg, res, iw, w, mem);
  }

  void Function::setup(const double** arg, double** res, int* iw, double* w,
                          int mem) const {
    (*this)->_setup(arg, res, iw, w, mem);
  }

  bool Function::spCanEvaluate(bool fwd) {
    if (fwd) {
      return (*this)->has_spfwd();
    } else {
      return (*this)->has_sprev();
    }
  }

  Function Function::derivative(int nfwd, int nadj) {
    // Quick return
    if (nfwd==0 && nadj==0) return *this;

    // Call self
    vector<MX> arg = mx_in();
    vector<MX> res = (*this)(arg);
    vector<MX> ret_in(arg), ret_out(res);

    // Number inputs and outputs
    int num_in = n_in();
    int num_out = n_out();

    // Forward sensitivities
    if (nfwd>0) {
      Function dfcn = forward(nfwd);
      arg = dfcn.mx_in();
      copy(ret_in.begin(), ret_in.begin()+num_in, arg.begin());
      copy(ret_out.begin(), ret_out.begin()+num_out, arg.begin()+num_in);
      ret_in.insert(ret_in.end(), arg.begin()+num_in+num_out, arg.end());
      res = dfcn(arg);
      vector<MX>::iterator it=res.begin();
      for (int d=0; d<nfwd; ++d)
        for (int i=0; i<num_out; ++i, ++it)
          *it = project(*it, sparsity_out(i));
      ret_out.insert(ret_out.end(), res.begin(), res.end());
    }

    // Adjoint sensitivities
    if (nadj>0) {
      Function dfcn = reverse(nadj);
      arg = dfcn.mx_in();
      copy(ret_in.begin(), ret_in.begin()+num_in, arg.begin());
      copy(ret_out.begin(), ret_out.begin()+num_out, arg.begin()+num_in);
      ret_in.insert(ret_in.end(), arg.begin()+num_in+num_out, arg.end());
      res = dfcn(arg);
      vector<MX>::iterator it=res.begin();
      for (int d=0; d<nadj; ++d)
        for (int i=0; i<num_in; ++i, ++it)
          *it = project(*it, sparsity_in(i));
      ret_out.insert(ret_out.end(), res.begin(), res.end());
    }

    // Name of return function
    stringstream ss;
    ss << "derivative_" << name() << "_" << nfwd << "_" << nadj;

    // Names of inputs
    vector<string> i_names;
    i_names.reserve(n_in()*(1+nfwd)+n_out()*nadj);
    const vector<string>& ischeme=(*this)->ischeme_;
    const vector<string>& oscheme=(*this)->oscheme_;

    // Nondifferentiated inputs
    for (int i=0; i<n_in(); ++i) {
      i_names.push_back("der_" + ischeme.at(i));
    }

    // Forward seeds
    for (int d=0; d<nfwd; ++d) {
      for (int i=0; i<n_in(); ++i) {
        ss.str(string());
        ss << "fwd" << d << "_" << ischeme.at(i);
        i_names.push_back(ss.str());
      }
    }

    // Adjoint seeds
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<n_out(); ++i) {
        ss.str(string());
        ss << "adj" << d << "_" << oscheme.at(i);
        i_names.push_back(ss.str());
      }
    }

    // Names of outputs
    vector<string> o_names;
    o_names.reserve(n_out()*(1+nfwd)+n_in()*nadj);

    // Nondifferentiated inputs
    for (int i=0; i<n_out(); ++i) {
      o_names.push_back("der_" + oscheme.at(i));
    }

    // Forward sensitivities
    for (int d=0; d<nfwd; ++d) {
      for (int i=0; i<n_out(); ++i) {
        ss.str(string());
        ss << "fwd" << d << "_" << oscheme.at(i);
        o_names.push_back(ss.str());
      }
    }

    // Adjoint sensitivities
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<n_in(); ++i) {
        ss.str(string());
        ss << "adj" << d << "_" << ischeme.at(i);
        o_names.push_back(ss.str());
      }
    }

    // Construct return function
    Function ret(ss.str(), ret_in, ret_out,
                 {{"input_scheme", i_names}, {"output_scheme", o_names}});

    // Consistency check for inputs
    int ind=0;
    for (int d=-1; d<nfwd; ++d) {
      for (int i=0; i<n_in(); ++i, ++ind) {
        if (ret.nnz_in(ind)!=0 && ret.sparsity_in(ind)!=sparsity_in(i)) {
          casadi_error("Incorrect sparsity for " << ret << " input " << ind << " \""
                       << i_names.at(ind) << "\". Expected " << size_in(i)
                       << " but got " << ret.size_in(ind));
        }
      }
    }
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<n_out(); ++i, ++ind) {
        if (ret.nnz_in(ind)!=0 && ret.sparsity_in(ind)!=sparsity_out(i)) {
          casadi_error("Incorrect sparsity for " << ret << " input " << ind <<
                       " \"" << i_names.at(ind) << "\". Expected " << size_out(i)
                       << " but got " << ret.size_in(ind));
        }
      }
    }

    // Consistency check for outputs
    ind=0;
    for (int d=-1; d<nfwd; ++d) {
      for (int i=0; i<n_out(); ++i, ++ind) {
        if (ret.nnz_out(ind)!=0 && ret.sparsity_out(ind)!=sparsity_out(i)) {
          casadi_error("Incorrect sparsity for " << ret << " output " << ind <<
                       " \"" <<  o_names.at(ind) << "\". Expected " << size_out(i)
                       << " but got " << ret.size_out(ind));
        }
      }
    }
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<n_in(); ++i, ++ind) {
        if (ret.nnz_out(ind)!=0 && ret.sparsity_out(ind)!=sparsity_in(i)) {
          casadi_error("Incorrect sparsity for " << ret << " output " << ind << " \""
                       << o_names.at(ind) << "\". Expected " << size_in(i)
                       << " but got " << ret.size_out(ind));
        }
      }
    }
    return ret;
  }

  Function Function::forward(int nfwd) {
    return (*this)->forward_old(nfwd);
  }

  Function Function::forward_new(int nfwd) {
    return (*this)->forward(nfwd);
  }

  Function Function::reverse(int nadj) {
    return (*this)->reverse_old(nadj);
  }

  Function Function::reverse_new(int nadj) {
    return (*this)->reverse(nadj);
  }

  void Function::set_forward(const Function& fcn, int nfwd) {
    (*this)->set_forward(fcn, nfwd);
  }

  void Function::set_reverse(const Function& fcn, int nadj) {
    (*this)->set_reverse(fcn, nadj);
  }

  void Function::print_dimensions(ostream &stream) const {
    (*this)->print_dimensions(stream);
  }

  void Function::print_options(ostream &stream) const {
    (*this)->print_options(stream);
  }

  void Function::print_option(const std::string &name, std::ostream &stream) const {
    (*this)->print_option(name, stream);
  }

  void Function::print_free(std::ostream &stream) const {
    (*this)->print_free(stream);
  }

  std::string Function::generate(const Dict& opts) {
    return generate(name(), opts);
  }

  std::string Function::generate(const string& fname, const Dict& opts) {
    CodeGenerator gen(fname, opts);
    gen.add(*this);
    return gen.generate();
  }

  std::string Function::generate_dependencies(const string& fname, const Dict& opts) {
    return (*this)->generate_dependencies(fname, opts);
  }

  void Function::checkInputs() const {
    return (*this)->checkInputs();
  }

  void Function::derivative(const DMVector& arg, DMVector& res,
                          const DMVectorVector& fseed, DMVectorVector& fsens,
                          const DMVectorVector& aseed, DMVectorVector& asens,
                          bool always_inline, bool never_inline) {
    call(arg, res, always_inline, never_inline);
    forward(arg, res, fseed, fsens, always_inline, never_inline);
    reverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  void Function::derivative(const SXVector& arg, SXVector& res,
                          const SXVectorVector& fseed, SXVectorVector& fsens,
                          const SXVectorVector& aseed, SXVectorVector& asens,
                          bool always_inline, bool never_inline) {
    call(arg, res, always_inline, never_inline);
    forward(arg, res, fseed, fsens, always_inline, never_inline);
    reverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  void Function::derivative(const MXVector& arg, MXVector& res,
                          const MXVectorVector& fseed, MXVectorVector& fsens,
                          const MXVectorVector& aseed, MXVectorVector& asens,
                          bool always_inline, bool never_inline) {
    call(arg, res, always_inline, never_inline);
    forward(arg, res, fseed, fsens, always_inline, never_inline);
    reverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  string Function::name() const {
    if (is_null()) {
      return "null";
    } else {
      return (*this)->name();
    }
  }

  bool Function::check_name(const std::string& name) {
    // Check if empty
    if (name.empty()) return false;

    // Check if keyword
    for (const char* kw : {"null", "jac", "hess"}) {
      if (name.compare(kw)==0) return false;
    }

    // Make sure that the first character is a letter
    auto it=name.begin();
    if (!std::isalpha(*it++)) return false;

    // Check remain_ing characters
    for (; it!=name.end(); ++it) {
      if (*it=='_') {
        // Make sure that the next character isn't also an underscore
        if (it+1!=name.end() && *(it+1)=='_') return false;
      } else {
        // Make sure alphanumeric
        if (!std::isalnum(*it)) return false;
      }
    }

    // Valid function name if reached this point
    return true;
  }

  string Function::fix_name(const string& name) {
    // Quick return if already valid name
    if (check_name(name)) return name;

    // If empty, name it "unnamed"
    if (name.empty()) return "unnamed";

    // Construct a sane name
    stringstream ss;

    // If the first character isn't a character, prepend an "a"
    if (!std::isalpha(name.front())) ss << "a";

    // Treat other characters
    bool previous_is_underscore = false;
    for (char c : name) {
      if (std::isalnum(c)) {
        // Alphanumeric characters
        ss << c;
        previous_is_underscore = false;
      } else if (!previous_is_underscore) {
        // Everything else becomes an underscore
        ss << '_';
        previous_is_underscore = true;
      }
    }

    // If name became a keyword, append 1
    for (const char* kw : {"null", "jac", "hess"}) {
      if (ss.str().compare(kw)==0) ss << "1";
    }

    return ss.str();
  }

  void Function::forward(const vector<MX>& arg, const vector<MX>& res,
                         const vector<vector<MX> >& fseed,
                         vector<vector<MX> >& fsens,
                         bool always_inline, bool never_inline) {
    (*this)->forward(arg, res, fseed, fsens, always_inline, never_inline);
  }

  void Function::reverse(const vector<MX>& arg, const vector<MX>& res,
                         const vector<vector<MX> >& aseed,
                         vector<vector<MX> >& asens,
                         bool always_inline, bool never_inline) {
    (*this)->reverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  void Function::forward(const vector<SX>& arg, const vector<SX>& res,
                         const vector<vector<SX> >& fseed,
                         vector<vector<SX> >& fsens,
                         bool always_inline, bool never_inline) {
    (*this)->forward(arg, res, fseed, fsens, always_inline, never_inline);
  }

  void Function::reverse(const vector<SX>& arg, const vector<SX>& res,
                         const vector<vector<SX> >& aseed,
                         vector<vector<SX> >& asens,
                         bool always_inline, bool never_inline) {
    (*this)->reverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  void Function::forward(const vector<DM>& arg, const vector<DM>& res,
                         const vector<vector<DM> >& fseed,
                         vector<vector<DM> >& fsens,
                         bool always_inline, bool never_inline) {
    (*this)->forward(arg, res, fseed, fsens, always_inline, never_inline);
  }

  void Function::reverse(const vector<DM>& arg, const vector<DM>& res,
                         const vector<vector<DM> >& aseed,
                         vector<vector<DM> >& asens,
                         bool always_inline, bool never_inline) {
    (*this)->reverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  vector<DM> Function::operator()(const vector<DM>& arg) {
    vector<DM> res;
    call(arg, res);
    return res;
  }

  vector<SX> Function::operator()(const vector<SX>& arg) {
    vector<SX> res;
    call(arg, res);
    return res;
  }

  vector<MX> Function::operator()(const vector<MX>& arg) {
    vector<MX> res;
    call(arg, res);
    return res;
  }

  template<typename M>
  void Function::_call(const std::map<string, M>& arg, std::map<string, M>& res,
                       bool always_inline, bool never_inline) {
    // Get default inputs
    vector<M> arg_v(n_in());
    for (int i=0; i<arg_v.size(); ++i) {
      arg_v[i] = default_in(i);
    }

    // Assign provided inputs
    for (auto&& e : arg) {
      arg_v.at(index_in(e.first)) = e.second;
    }

    // Make call
    vector<M> res_v;
    call(arg_v, res_v, always_inline, never_inline);

    // Save to map
    res.clear();
    for (int i=0; i<res_v.size(); ++i) {
      res[name_out(i)] = res_v[i];
    }
  }

  const DMDict Function::operator()(const DMDict& arg) {
    DMDict res;
    call(arg, res);
    return res;
  }

  const SXDict Function::operator()(const SXDict& arg) {
    SXDict res;
    call(arg, res);
    return res;
  }

  const MXDict Function::operator()(const MXDict& arg) {
    MXDict res;
    call(arg, res);
    return res;
  }

  void Function::call(const DMDict& arg, DMDict& res,
                      bool always_inline, bool never_inline) {
    return _call(arg, res, always_inline, never_inline);
  }

  void Function::call(const SXDict& arg, SXDict& res,
                      bool always_inline, bool never_inline) {
    return _call(arg, res, always_inline, never_inline);
  }

  void Function::call(const MXDict& arg, MXDict& res,
                      bool always_inline, bool never_inline) {
    return _call(arg, res, always_inline, never_inline);
  }

  double Function::default_in(int ind) const {
    return (*this)->default_in(ind);
  }

  void Function::operator()(const double** arg, double** res, int* iw, double* w, int mem) const {
    (*this)->_eval(arg, res, iw, w, mem);
  }

  void Function::operator()(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) const {
    (*this)->eval_sx(arg, res, iw, w, mem);
  }

  const SX Function::sx_in(int ind) const {
    return (*this)->sx_in(ind);
  }

  const SX Function::sx_out(int ind) const {
    return (*this)->sx_out(ind);
  }

  const vector<SX> Function::sx_in() const {
    return (*this)->sx_in();
  }

  const vector<SX> Function::sx_out() const {
    return (*this)->sx_out();
  }

  const MX Function::mx_in(int ind) const {
    return (*this)->mx_in(ind);
  }

  const MX Function::mx_out(int ind) const {
    return (*this)->mx_out(ind);
  }

  const vector<MX> Function::mx_in() const {
    return (*this)->mx_in();
  }

  const vector<MX> Function::mx_out() const {
    return (*this)->mx_out();
  }

  string Function::type_name() const {
    return (*this)->type_name();
  }

  bool Function::is_a(const string& type, bool recursive) const {
    return (*this)->is_a(type, recursive);
  }

  vector<SX> Function::free_sx() const {
    return (*this)->free_sx();
  }

  vector<MX> Function::free_mx() const {
    return (*this)->free_mx();
  }

  bool Function::has_free() const {
    return (*this)->has_free();
  }

  void Function::generate_lifted(Function& vdef_fcn, Function& vinit_fcn) {
    (*this)->generate_lifted(vdef_fcn, vinit_fcn);
  }

  int Function::getAlgorithmSize() const {
    return (*this)->getAlgorithmSize();
  }

  int Function::getWorkSize() const {
    return (*this)->getWorkSize();
  }

  int Function::getAtomicOperation(int k) const {
    return (*this)->getAtomicOperation(k);
  }

  pair<int, int> Function::getAtomicInput(int k) const {
    return (*this)->getAtomicInput(k);
  }

  double Function::getAtomicInputReal(int k) const {
    return (*this)->getAtomicInputReal(k);
  }

  int Function::getAtomicOutput(int k) const {
    return (*this)->getAtomicOutput(k);
  }

  int Function::n_nodes() const {
    return (*this)->n_nodes();
  }

  int Function::checkout() {
    return (*this)->checkout();
  }

  void Function::release(int mem) {
    (*this)->release(mem);
  }

  void* Function::memory(int ind) const {
    return (*this)->memory(ind);
  }

  void Function::assert_size_in(int i, int nrow, int ncol) const {
    casadi_assert_message(size1_in(i)==nrow && size2_in(i)==ncol,
                          "Incorrect shape for " << *this << " input " << i << " \""
                          << name_in(i) << "\". Expected " << nrow << "-by-" << ncol
                          << " but got " << size1_in(i) <<  "-by-" << size2_in(i));

  }

  void Function::assert_size_out(int i, int nrow, int ncol) const {
    casadi_assert_message(size1_out(i)==nrow && size2_out(i)==ncol,
                          "Incorrect shape for " << *this << " output " << i << " \""
                          << name_out(i) << "\". Expected " << nrow << "-by-" << ncol
                          << " but got " << size1_out(i) <<  "-by-" << size2_out(i));
  }

  Function Function::
  factory(const std::string& name,
          const std::vector<std::string>& s_in,
          const std::vector<std::string>& s_out,
          const AuxOut& aux,
          const Dict& opts) const {
     return (*this)->factory(name, s_in, s_out, aux, opts);
  }

  vector<bool> Function::nl_var(const string& s_in,
                                const vector<string>& s_out) const {
    return (*this)->nl_var(s_in, s_out);
  }

  std::vector<std::string> Function::get_function() const {
    return (*this)->get_function();
  }

  Function Function::get_function(const std::string &name) const {
    return (*this)->get_function(name);
  }

  bool Function::has_function(const std::string& fname) const {
    return (*this)->has_function(fname);
  }

  Function Function::oracle() const {
    return (*this)->oracle();
  }

} // namespace casadi
