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
#include "external.hpp"
#include "switch.hpp"
#include "kernel_sum.hpp"
#include "jit.hpp"
#include "nlpsol.hpp"
#include "qpsol.hpp"

#include <typeinfo>

using namespace std;

namespace casadi {

  Function::Function() {
  }

  Function::~Function() {
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

  Function Function::mapaccum(const string& name, int N, const Dict& opts) const {
    vector<bool> accum_input(n_in(), false);
    accum_input[0] = true;
    vector<int> accum_output(1, 0);
    return mapaccum(name, N, accum_input, accum_output, false, opts);
  }

  Function Function::mapaccum(const string& name, int n,
                              const vector<bool>& input_accum,
                              const vector<int>& output_accum,
                              bool reverse,
                              const Dict& opts) const {
    Function ret;
    ret.assignNode(new Mapaccum(name, *this, n, input_accum, output_accum, reverse));
    ret->construct(opts);
    return ret;
  }

  Function Function::map(const string& name, int n, const Dict& opts) const {
    Function ret;
    ret.assignNode(MapBase::create(name, *this, n, opts));
    ret->construct(opts);
    return ret;
  }

  Function Function::map(const string& name,
                         int n,
                         const vector<bool> &repeat_in,
                         const vector<bool> &repeat_out,
                         const Dict& opts) const {

    Function ret;
    ret.assignNode(new MapReduce(name, *this, n, repeat_in, repeat_out));
    ret->construct(opts);
    return ret;
  }

  vector<vector<MX> > Function::map(const vector<vector<MX> > &x,
                                    const string& parallelization) {
    return (*this)->map_mx(x, parallelization);
  }

  vector<MX> Function::map(const vector< MX > &x,
                           const string& parallelization) {
    return (*this)->map_mx(x, parallelization);
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

  void Function::addMonitor(const string& mon) {
    (*this)->monitors_.insert(mon);
  }

  void Function::removeMonitor(const string& mon) {
    (*this)->monitors_.erase(mon);
  }

  Dict Function::stats(int mem) const {
    return (*this)->mem_.at(mem)->get_stats();
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

  string Function::description_in(int ind) const {
    return (*this)->description_in(ind);
  }

  string Function::description_out(int ind) const {
    return (*this)->description_out(ind);
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
    (*const_cast<Function*>(this))->spFwd(arg, res, iw, w, mem);
  }

  void Function::rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    (*this)->spAdj(arg, res, iw, w, mem);
  }

  void Function::set_work(const double**& arg, double**& res, int*& iw, double*& w,
                          int mem) const {
    (*this)->set_work(*(*this)->mem_.at(mem), arg, res, iw, w);
  }

  void Function::set_temp(const double** arg, double** res, int* iw, double* w,
                          int mem) const {
    (*this)->set_temp(*(*this)->mem_.at(mem), arg, res, iw, w);
  }

  void Function::setup(const double** arg, double** res, int* iw, double* w,
                          int mem) const {
    (*this)->setup(*(*this)->mem_.at(mem), arg, res, iw, w);
  }

  bool Function::spCanEvaluate(bool fwd) {
    return (*this)->spCanEvaluate(fwd);
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
    return (*this)->forward(nfwd);
  }

  Function Function::reverse(int nadj) {
    return (*this)->reverse(nadj);
  }

  void Function::set_forward(const Function& fcn, int nfwd) {
    (*this)->set_forward(fcn, nfwd);
  }

  void Function::set_reverse(const Function& fcn, int nadj) {
    (*this)->set_reverse(fcn, nadj);
  }

  void Function::printDimensions(ostream &stream) const {
    (*this)->printDimensions(stream);
  }

  void Function::printOptions(ostream &stream) const {
    (*this)->printOptions(stream);
  }

  void Function::printOption(const std::string &name, std::ostream &stream) const {
    (*this)->printOption(name, stream);
  }

  void Function::generate(const Dict& opts) {
    generate(getSanitizedName(), opts);
  }

  void Function::generate(const string& fname, const Dict& opts) {
    CodeGenerator gen(opts);
    gen.add(*this, fname);
    gen.generate(fname);
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
      return "NULL";
    } else {
      return (*this)->name();
    }
  }

  string Function::getSanitizedName() const {
    return (*this)->getSanitizedName();
  }

  string Function::sanitizeName(const string& name) {
    return FunctionInternal::sanitizeName(name);
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

  vector<DM> Function::operator()(const vector<DM>& arg,
                                  bool always_inline, bool never_inline) {
    vector<DM> res;
    call(arg, res, always_inline, never_inline);
    return res;
  }

  vector<SX> Function::operator()(const vector<SX>& arg,
                                  bool always_inline, bool never_inline) {
    vector<SX> res;
    call(arg, res, always_inline, never_inline);
    return res;
  }

  vector<MX> Function::operator()(const vector<MX>& arg,
                                  bool always_inline, bool never_inline) {
    vector<MX> res;
    call(arg, res, always_inline, never_inline);
    return res;
  }

  template<typename M>
  const std::map<string, M>
  Function::callMap(const std::map<string, M>& arg, bool always_inline, bool never_inline) {
    // Get default inputs
    vector<M> v(n_in());
    for (int i=0; i<v.size(); ++i) {
      v[i] = default_in(i);
    }

    // Assign provided inputs
    for (typename std::map<string, M>::const_iterator i=arg.begin(); i!=arg.end(); ++i) {
      v.at(index_in(i->first)) = i->second;
    }

    // Make call
    v = (*this)(v, always_inline, never_inline);

    // Save to map
    std::map<string, M> ret;
    for (int i=0; i<v.size(); ++i) {
      ret[name_out(i)] = v[i];
    }
    return ret;
  }

  const DMDict Function::operator()(const DMDict& arg, bool always_inline,
                                    bool never_inline) {
    return callMap(arg, always_inline, never_inline);
  }

  const SXDict Function::operator()(const SXDict& arg, bool always_inline,
                                    bool never_inline) {
    return callMap(arg, always_inline, never_inline);
  }

  const MXDict Function::operator()(const MXDict& arg, bool always_inline,
                                    bool never_inline) {
    return callMap(arg, always_inline, never_inline);
  }

  double Function::default_in(int ind) const {
    return (*this)->default_in(ind);
  }

  void Function::operator()(const double** arg, double** res, int* iw, double* w, int mem) const {
    (*this)->eval(arg, res, iw, w, mem);
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

  SX Function::free_sx() const {
    return (*this)->free_sx();
  }

  vector<MX> Function::free_mx() const {
    return (*this)->free_mx();
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

  int Function::countNodes() const {
    return (*this)->countNodes();
  }

  Function external(const string& name, const Dict& opts) {
    Function ret;
    ret.assignNode(External::create("./" + name + ".so", name));
    ret->construct(opts);
    return ret;
  }

  Function external(const string& name, const string& bin_name,
                              const Dict& opts) {
    Function ret;
    ret.assignNode(External::create(bin_name, name));
    ret->construct(opts);
    return ret;
  }

  Function external(const string& name, const Compiler& compiler,
                              const Dict& opts) {
    Function ret;
    ret.assignNode(External::create(compiler, name));
    ret->construct(opts);
    return ret;
  }

  XProblem::XProblem(const SXProblem& d) : sx_p(new SXProblem(d)), is_sx(true) {
  }

  XProblem::XProblem(const MXProblem& d) : mx_p(new MXProblem(d)), is_sx(false) {
  }

  XProblem::~XProblem() {
    if (is_sx) {
      delete sx_p;
    } else {
      delete mx_p;
    }
  }

  XProblem::XProblem(const XProblem& d) : is_sx(d.is_sx) {
    if (d.is_sx) {
      sx_p = new SXProblem(*d.sx_p);
    } else {
      mx_p = new MXProblem(*d.mx_p);
    }
  }

  XProblem& XProblem::operator=(const XProblem& d) {
    if (&d!=this) {
      // Delete the previous object
      if (is_sx) {
        delete sx_p;
      } else {
        delete mx_p;
      }
      // Assign
      is_sx = d.is_sx;
      if (is_sx) {
        sx_p = new SXProblem(*d.sx_p);
      } else {
        mx_p = new MXProblem(*d.mx_p);
      }
    }
    return *this;
  }

  XProblem::operator const SXProblem&() const {
    casadi_assert(is_sx);
    return *sx_p;
  }

  XProblem::operator const MXProblem&() const {
    casadi_assert(!is_sx);
    return *mx_p;
  }

  Function jit(const std::string& name, int n_in, int n_out,
               const std::string& body, const Dict& opts) {
    Function ret;
    ret.assignNode(new Jit(name, n_in, n_out, body, opts));
    ret->construct(opts);
    return ret;
  }

} // namespace casadi

