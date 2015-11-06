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
#include "ivpsol.hpp"
#include "qpsol.hpp"
#include "nlpsol.hpp"
#include "nlsol.hpp"
#include "linsol.hpp"


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
    setOption(opts);
    init();
  }

  void Function::construct(const string& name,
                           const vector<MX>& arg, const vector<MX>& res,
                           const Dict& opts) {
    assignNode(new MXFunction(name, arg, res));
    setOption(opts);
    init();
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

  const FunctionInternal* Function::operator->() const {
    return static_cast<const FunctionInternal*>(OptionsFunctionality::operator->());
  }

  FunctionInternal* Function::operator->() {
    return static_cast<FunctionInternal*>(OptionsFunctionality::operator->());
  }

  void Function::call(const vector<DMatrix> &arg, vector<DMatrix> &res,
                      bool always_inline, bool never_inline) {
    if (!matchingArg(arg))
      return call(replaceArg(arg), res, always_inline, never_inline);
    (*this)->call(arg, res, always_inline, never_inline);
  }

  void Function::call(const vector<SX> &arg, vector<SX>& res,
                      bool always_inline, bool never_inline) {
    if (!matchingArg(arg))
      return call(replaceArg(arg), res, always_inline, never_inline);
    (*this)->call(arg, res, always_inline, never_inline);
  }

  void Function::call(const vector<MX> &arg, vector<MX>& res,
                      bool always_inline, bool never_inline) {
    if (!matchingArg(arg))
      return call(replaceArg(arg), res, always_inline, never_inline);
    (*this)->call(arg, res, always_inline, never_inline);
  }

  vector<const double*> Function::buf_in(Function::VecArg arg) const {
    casadi_assert(arg.size()==n_in());
    auto arg_it=arg.begin();
    vector<const double*> buf_arg(sz_arg());
    for (unsigned int i=0; i<arg.size(); ++i) {
      casadi_assert(arg_it->size()==nnz_in(i));
      buf_arg[i] = getPtr(*arg_it++);
    }
    return buf_arg;
  }

  vector<const double*> Function::buf_in(Function::L1dArg arg) const {
    casadi_assert(arg.size()==n_in());
    auto arg_it=arg.begin();
    vector<const double*> buf_arg(sz_arg());
    for (unsigned int i=0; i<arg.size(); ++i) {
      casadi_assert(arg_it->size()==nnz_in(i));
      buf_arg[i] = getPtr(*arg_it++);
    }
    return buf_arg;
  }

  vector<double*> Function::buf_out(Function::VecRes res) const {
    res.resize(n_out());
    auto res_it=res.begin();
    vector<double*> buf_res(sz_res());
    for (unsigned int i=0; i<res.size(); ++i) {
      res_it->resize(nnz_out(i));
      buf_res[i] = getPtr(*res_it++);
    }
    return buf_res;
  }

  vector<double*> Function::buf_out(Function::L1dRes res) const {
    casadi_assert(res.size()==n_out());
    auto res_it=res.begin();
    vector<double*> buf_res(sz_res());
    for (unsigned int i=0; i<res.size(); ++i) {
      casadi_assert(*res_it!=0);
      (*res_it)->resize(nnz_out(i));
      buf_res[i] = getPtr(**res_it++);
    }
    return buf_res;
  }

  vector<const double*> Function::buf_in(Function::MapArg arg) const {
    // Return value (RVO)
    vector<const double*> ret(sz_arg());

    // Get default values
    for (int i=0; i<n_in(); ++i) ret[i] = &default_in(i);

    // Read inputs
    for (auto i=arg.begin(); i!=arg.end(); ++i) {
      int ind = index_in(i->first);
      casadi_assert(i->second.size()==nnz_in(ind));
      ret[ind] = getPtr(i->second);
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
      ret[ind] = getPtr(i->second);
    }

    return ret;
  }

  vector<const double*> Function::buf_in(Function::L2dArg arg) const {
    // Return value (RVO)
    vector<const double*> ret(sz_arg());

    // Get default values
    for (int i=0; i<n_in(); ++i) ret[i] = &default_in(i);

    // Read inputs
    for (auto i=arg.begin(); i!=arg.end(); ++i) {
      int ind = index_in(i->first);
      casadi_assert(i->second.size()==nnz_in(ind));
      ret[ind] = getPtr(i->second);
    }

    return ret;
  }

  vector<double*> Function::buf_out(Function::L2dRes res) const {
    // Return value (RVO)
    vector<double*> ret(sz_res(), 0);

    // Read outputs
    for (auto i=res.begin(); i!=res.end(); ++i) {
      int ind = index_out(i->first);
      casadi_assert(i->second!=0);
      i->second->resize(nnz_out(ind));
      ret[ind] = getPtr(*i->second);
    }

    return ret;
  }

  void Function::operator()(vector<const double*> arg, vector<double*> res) {
    casadi_assert(arg.size()>=sz_arg());
    casadi_assert(res.size()>=sz_res());
    vector<int> buf_iw(sz_iw());
    vector<double> buf_w(sz_w());

    // Evaluate memoryless
    (*this)(0, getPtr(arg), getPtr(res), getPtr(buf_iw), getPtr(buf_w));
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
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  Function Function::map(const string& name, int n, const Dict& opts) const {
    Function ret;
    ret.assignNode(MapBase::create(name, *this, n, opts));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  Function Function::map(const string& name,
                         int n,
                         const vector<bool> &repeat_in,
                         const vector<bool> &repeat_out,
                         const Dict& opts) const {

    Function ret;
    ret.assignNode(new MapReduce(name, *this, n, repeat_in, repeat_out));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  vector<vector<MX> > Function::map(const vector<vector<MX> > &x,
                                    const string& parallelization) {
    if (x.empty()) return x;

    // Check if arguments match
    int n = x.size();
    bool matching=true;
    for (int i=0; i<n; ++i) {
      matching = matchingArg(x[i]) && matching; // non-short circuiting
    }

    // Replace arguments if needed
    if (!matching) {
      vector<vector<MX> > x_new(n);
      for (int i=0; i<n; ++i) x_new[i] = replaceArg(x[i]);
      return map(x_new, parallelization);
    }

    vector< vector<MX> > trans = swapIndices(x);
    vector< MX > x_cat(trans.size());
    for (int i=0;i<trans.size();++i) {
      x_cat[i] = horzcat(trans[i]);
    }

    // Call the internal function
    vector< MX > ret_cat = map(x_cat, parallelization);
    vector< vector<MX> > ret;

    for (int i=0;i<ret_cat.size();++i) {
      ret.push_back(horzsplit(ret_cat[i], output(i).size2()));
    }
    return swapIndices(ret);

  }

  vector<MX> Function::map(const vector< MX > &x,
                                    const string& parallelization) {
    if (x.empty()) return x;

    // Replace arguments if needed
    if (!matchingArg(x, true)) {
      vector< MX > x_new = replaceArg(x, true);
      return map(x_new, parallelization);
    }

    int n = 1;
    for (int i=0;i<x.size();++i) {
      n = max(x[i].size2()/input(i).size2(), n);
    }

    vector<bool> repeat_n;
    for (int i=0;i<x.size();++i) {
      repeat_n.push_back(x[i].size2()/input(i).size2()==n);
    }

    bool repeated = true;
    for (int i=0;i<repeat_n.size();++i) {
      repeated &= repeat_n[i];
    }

    // Call the internal function

    Dict options = {{"parallelization", parallelization}};

    Function ms;
    if (repeated) {
      ms = map("map", n, options);
    } else {
      ms = map("mapsum", n, repeat_n, vector<bool>(n_out(), true), options);
    }
    // Call the internal function
    return ms(x);
  }

  vector<MX> Function::mapsum(const vector< MX > &x,
                                    const string& parallelization) {
    if (x.empty()) return x;

    // Replace arguments if needed
    if (!matchingArg(x, true)) {
      vector< MX > x_new = replaceArg(x, true);
      return mapsum(x_new, parallelization);
    }

    int n = 1;
    for (int i=0;i<x.size();++i) {
      n = max(x[i].size2()/input(i).size2(), n);
    }

    vector<bool> repeat_n;
    for (int i=0;i<x.size();++i) {
      repeat_n.push_back(x[i].size2()/input(i).size2()==n);
    }

    Dict options = {{"parallelization", parallelization}};
    Function ms = map("mapsum", n, repeat_n,
                      vector<bool>(n_out(), false), options);

    // Call the internal function
    return ms(x);
  }

  Function Function::conditional(const string& name, const vector<Function>& f,
                                 const Function& f_def, const Dict& opts) {
    Function ret;
    ret.assignNode(new Switch(name, f, f_def));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  Function Function::if_else(const string& name, const Function& f_true,
                             const Function& f_false, const Dict& opts) {
    Function ret;
    ret.assignNode(new Switch(name, vector<Function>(1, f_false), f_true));
    ret.setOption(opts);
    ret.init();
    return ret;
  }


  Function Function::kernel_sum(const string& name,
                                const pair<int, int> & size,
                                double r, int n,
                                const Dict& opts) const {
    Function ret;
    ret.assignNode(new KernelSum(name, *this, size, r, n));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  void Function::evaluate() {
    (*this)->evaluate();
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

  const Dict & Function::getStats() const {
    return (*this)->getStats();
  }

  GenericType Function::getStat(const string& name) const {
    return (*this)->getStat(name);
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

  const Matrix<double>& Function::input(int i) const {
    return (*this)->input(i);
  }

  const Matrix<double>& Function::input(const string &iname) const {
    return (*this)->input(iname);
  }

  Sparsity Function::sparsity_in(int ind) const {
    return (*this)->sparsity_in(ind);
  }

  Sparsity Function::sparsity_in(const string &iname) const {
    return (*this)->sparsity_in(iname);
  }

  Sparsity Function::sparsity_out(int ind) const {
    return (*this)->sparsity_out(ind);
  }

  Sparsity Function::sparsity_out(const string &iname) const {
    return (*this)->sparsity_out(iname);
  }

  Matrix<double>& Function::input(int i) {
    return (*this)->input(i);
  }

  Matrix<double>& Function::input(const string &iname) {
    return (*this)->input(iname);
  }

  const Matrix<double>& Function::output(int i) const {
    return (*this)->output(i);
  }

  const Matrix<double>& Function::output(const string &oname) const {
    return (*this)->output(oname);
  }

  Matrix<double>& Function::output(int i) {
    return (*this)->output(i);
  }

  Matrix<double>& Function::output(const string &oname) {
    return (*this)->output(oname);
  }

  void Function::spEvaluate(bool fwd) {
    (*this)->spEvaluate(fwd);
  }

  void Function::sz_work(size_t& sz_arg, size_t& sz_res, size_t& sz_iw, size_t& sz_w) const {
    (*this)->sz_work(sz_arg, sz_res, sz_iw, sz_w);
  }

  size_t Function::sz_arg() const { return (*this)->sz_arg();}

  size_t Function::sz_res() const { return (*this)->sz_res();}

  size_t Function::sz_iw() const { return (*this)->sz_w();}

  size_t Function::sz_w() const { return (*this)->sz_w();}

  void* Function::alloc_mem() { return (*this)->alloc_mem();}

  void Function::free_mem(void* mem) { (*this)->free_mem(mem);}

  MemBlock Function::alloc() const {
    casadi_assert(!isNull());
    const Function& f = *this;
    return MemBlock(f);
  }

  void Function::operator()(void* mem, const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    (*this)->spFwdSwitch(mem, arg, res, iw, w);
  }

  void Function::rev(void* mem, bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    (*this)->spAdjSwitch(mem, arg, res, iw, w);
  }

  bool Function::spCanEvaluate(bool fwd) {
    return (*this)->spCanEvaluate(fwd);
  }

  void Function::spInit(bool fwd) {
    (*this)->spInit(fwd);
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
          *it = project(*it, output(i).sparsity());
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
          *it = project(*it, input(i).sparsity());
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
        if (ret.input(ind).nnz()!=0 && ret.input(ind).sparsity()!=input(i).sparsity()) {
          casadi_error("Incorrect sparsity for " << ret << " input " << ind << " \""
                       << i_names.at(ind) << "\". Expected " << input(i).dim()
                       << " but got " << ret.input(ind).dim());
        }
      }
    }
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<n_out(); ++i, ++ind) {
        if (ret.input(ind).nnz()!=0 && ret.input(ind).sparsity()!=output(i).sparsity()) {
          casadi_error("Incorrect sparsity for " << ret << " input " << ind <<
                       " \"" << i_names.at(ind) << "\". Expected " << output(i).dim()
                       << " but got " << ret.input(ind).dim());
        }
      }
    }

    // Consistency check for outputs
    ind=0;
    for (int d=-1; d<nfwd; ++d) {
      for (int i=0; i<n_out(); ++i, ++ind) {
        if (ret.output(ind).nnz()!=0 && ret.output(ind).sparsity()!=output(i).sparsity()) {
          casadi_error("Incorrect sparsity for " << ret << " output " << ind <<
                       " \"" <<  o_names.at(ind) << "\". Expected " << output(i).dim()
                       << " but got " << ret.output(ind).dim());
        }
      }
    }
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<n_in(); ++i, ++ind) {
        if (ret.output(ind).nnz()!=0 && ret.output(ind).sparsity()!=input(i).sparsity()) {
          casadi_error("Incorrect sparsity for " << ret << " output " << ind << " \""
                       << o_names.at(ind) << "\". Expected " << input(i).dim()
                       << " but got " << ret.output(ind).dim());
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

  void Function::derivative(const DMatrixVector& arg, DMatrixVector& res,
                          const DMatrixVectorVector& fseed, DMatrixVectorVector& fsens,
                          const DMatrixVectorVector& aseed, DMatrixVectorVector& asens,
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
    return (*this)->name_;
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
    checkArg(arg);
    checkRes(res);
    if (!matchingFwdSeed(fseed)) {
      return forward(arg, res, replaceFwdSeed(fseed), fsens,
                     always_inline, never_inline);
    }
    (*this)->forward(arg, res, fseed, fsens, always_inline, never_inline);
  }

  void Function::reverse(const vector<MX>& arg, const vector<MX>& res,
                         const vector<vector<MX> >& aseed,
                         vector<vector<MX> >& asens,
                         bool always_inline, bool never_inline) {
    checkArg(arg);
    checkRes(res);
    if (!matchingAdjSeed(aseed)) {
      return reverse(arg, res, replaceAdjSeed(aseed), asens,
                     always_inline, never_inline);
    }
    (*this)->reverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  void Function::forward(const vector<SX>& arg, const vector<SX>& res,
                         const vector<vector<SX> >& fseed,
                         vector<vector<SX> >& fsens,
                         bool always_inline, bool never_inline) {
    checkArg(arg);
    checkRes(res);
    if (!matchingFwdSeed(fseed)) {
      return forward(arg, res, replaceFwdSeed(fseed), fsens,
                     always_inline, never_inline);
    }
    (*this)->forward(arg, res, fseed, fsens, always_inline, never_inline);
  }

  void Function::reverse(const vector<SX>& arg, const vector<SX>& res,
                         const vector<vector<SX> >& aseed,
                         vector<vector<SX> >& asens,
                         bool always_inline, bool never_inline) {
    checkArg(arg);
    checkRes(res);
    if (!matchingAdjSeed(aseed)) {
      return reverse(arg, res, replaceAdjSeed(aseed), asens,
                     always_inline, never_inline);
    }
    (*this)->reverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  void Function::forward(const vector<DMatrix>& arg, const vector<DMatrix>& res,
                         const vector<vector<DMatrix> >& fseed,
                         vector<vector<DMatrix> >& fsens,
                         bool always_inline, bool never_inline) {
    checkArg(arg);
    checkRes(res);
    if (!matchingFwdSeed(fseed)) {
      return forward(arg, res, replaceFwdSeed(fseed), fsens,
                     always_inline, never_inline);
    }
    (*this)->forward(arg, res, fseed, fsens, always_inline, never_inline);
  }

  void Function::reverse(const vector<DMatrix>& arg, const vector<DMatrix>& res,
                         const vector<vector<DMatrix> >& aseed,
                         vector<vector<DMatrix> >& asens,
                         bool always_inline, bool never_inline) {
    checkArg(arg);
    checkRes(res);
    if (!matchingAdjSeed(aseed)) {
      return reverse(arg, res, replaceAdjSeed(aseed), asens,
                     always_inline, never_inline);
    }
    (*this)->reverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  vector<DMatrix> Function::operator()(const vector<DMatrix>& arg,
                                            bool always_inline, bool never_inline) {
    vector<DMatrix> res;
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

  const DMatrixDict Function::operator()(const DMatrixDict& arg, bool always_inline,
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

  inline bool checkMat(const Sparsity& arg, const Sparsity& inp, bool hcat=false) {
    return arg.size()==inp.size() || arg.is_empty() || arg.is_scalar() ||
      (inp.size2()==arg.size1() && inp.size1()==arg.size2()
       && (arg.is_column() || inp.is_column())) ||
      (hcat && arg.size1()==inp.size1() && arg.size2() % inp.size2()==0);
  }

  template<typename M>
  void Function::checkArg(const vector<M>& arg, bool hcat) const {
    int n_in = this->n_in();
    casadi_assert_message(arg.size()==n_in, "Incorrect number of inputs: Expected "
                          << n_in << ", got " << arg.size());
    for (int i=0; i<n_in; ++i) {
      casadi_assert_message(checkMat(arg[i].sparsity(), input(i).sparsity(), hcat),
                            "Input " << i << " has mismatching shape. Expected "
                            << input(i).size() << ", got " << arg[i].size());
    }
  }

  template<typename M>
  void Function::checkRes(const vector<M>& res) const {
    int n_out = this->n_out();
    casadi_assert_message(res.size()==n_out, "Incorrect number of outputs: Expected "
                          << n_out << ", got " << res.size());
    for (int i=0; i<n_out; ++i) {
      casadi_assert_message(checkMat(res[i].sparsity(), output(i).sparsity()),
                            "Output " << i << " has mismatching shape. Expected "
                            << output(i).size() << ", got " << res[i].size());
    }
  }

  template<typename M>
  void Function::checkFwdSeed(const vector<vector<M> >& fseed) const {
    int n_in = this->n_in();
    for (int d=0; d<fseed.size(); ++d) {
      casadi_assert_message(fseed[d].size()==n_in,
                            "Incorrect number of forward seeds for direction " << d
                            << ": Expected " << n_in << ", got " << fseed[d].size());
      for (int i=0; i<n_in; ++i) {
        casadi_assert_message(checkMat(fseed[d][i].sparsity(), input(i).sparsity()),
                              "Forward seed " << i << " for direction " << d
                              << " has mismatching shape. Expected " << input(i).size()
                              << ", got " << fseed[d][i].size());
      }
    }
  }

  template<typename M>
  void Function::checkAdjSeed(const vector<vector<M> >& aseed) const {
    int n_out = this->n_out();
    for (int d=0; d<aseed.size(); ++d) {
      casadi_assert_message(aseed[d].size()==n_out,
                            "Incorrect number of adjoint seeds for direction " << d
                            << ": Expected " << n_out << ", got " << aseed[d].size());
      for (int i=0; i<n_out; ++i) {
        casadi_assert_message(checkMat(aseed[d][i].sparsity(), output(i).sparsity()),
                              "Adjoint seed " << i << " for direction " << d
                              << " has mismatching shape. Expected " << output(i).size()
                              << ", got " << aseed[d][i].size());
      }
    }
  }

  template<typename M>
  bool Function::matchingArg(const vector<M>& arg, bool hcat) const {
    checkArg(arg, hcat);
    int n_in = this->n_in();
    for (int i=0; i<n_in; ++i) {
      if (hcat) {
        if (arg.at(i).size1()!=input(i).size1()) return false;
        if (arg.at(i).size2() % input(i).size2()!=0 || arg.at(i).size2()==0) return false;
      } else {
        if (arg.at(i).size()!=input(i).size()) return false;
      }
    }
    return true;
  }

  template<typename M>
  bool Function::matchingRes(const vector<M>& res) const {
    checkRes(res);
    int n_out = this->n_out();
    for (int i=0; i<n_out; ++i) {
      if (res.at(i).size()!=output(i).size()) return false;
    }
    return true;
  }

  template<typename M>
  bool Function::matchingFwdSeed(const vector<vector<M> >& fseed) const {
    checkFwdSeed(fseed);
    for (int d=0; d<fseed.size(); ++d) {
      if (!matchingArg(fseed[d])) return false;
    }
    return true;
  }

  template<typename M>
  bool Function::matchingAdjSeed(const vector<vector<M> >& aseed) const {
    checkAdjSeed(aseed);
    for (int d=0; d<aseed.size(); ++d) {
      if (!matchingRes(aseed[d])) return false;
    }
    return true;
  }

  template<typename M>
  M replaceMat(const M& arg, const Sparsity& inp, bool hcat=false) {
    if (arg.size()==inp.size()) {
      // Matching dimensions already
      return arg;
    } else if (hcat && arg.size1()==inp.size1() && arg.size2() % inp.size2()==0
                    && arg.size2() >=0) {
      // Matching horzcat dimensions
      return arg;
    } else if (arg.is_empty()) {
      // Empty matrix means set zero
      return M(inp.size());
    } else if (arg.is_scalar()) {
      // Scalar assign means set all
      return M(inp, arg);
    } else {
      // Assign vector with transposing
      casadi_assert(arg.size1()==inp.size2() && arg.size2()==inp.size1()
                    && (arg.is_column() || inp.is_column()));
      return arg.T();
    }
  }

  template<typename M>
  vector<M> Function::replaceArg(const vector<M>& arg, bool hcat) const {
    vector<M> r(arg.size());
    for (int i=0; i<r.size(); ++i) r[i] = replaceMat(arg[i], input(i).sparsity(), hcat);
    return r;
  }

  template<typename M>
  vector<M> Function::replaceRes(const vector<M>& res) const {
    vector<M> r(res.size());
    for (int i=0; i<r.size(); ++i) r[i] = replaceMat(res[i], output(i).sparsity());
    return r;
  }

  template<typename M>
  vector<vector<M> >
  Function::replaceFwdSeed(const vector<vector<M> >& fseed) const {
    vector<vector<M> > r(fseed.size());
    for (int d=0; d<r.size(); ++d) r[d] = replaceArg(fseed[d]);
    return r;
  }

  template<typename M>
  vector<vector<M> >
  Function::replaceAdjSeed(const vector<vector<M> >& aseed) const {
    vector<vector<M> > r(aseed.size());
    for (int d=0; d<r.size(); ++d) r[d] = replaceRes(aseed[d]);
    return r;
  }

  const double& Function::default_in(int ind) const {
    return (*this)->default_in(ind);
  }

  void Function::operator()(void* mem, const double** arg, double** res, int* iw, double* w) {
    (*this)->eval(mem, arg, res, iw, w);
  }

  void Function::operator()(void* mem, const SXElem** arg, SXElem** res, int* iw, SXElem* w) {
    (*this)->evalSX(mem, arg, res, iw, w);
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

  void Function::init() {
    (*this)->init();
    (*this)->finalize();
  }

  Function Function::external(const string& name, const Dict& opts) {
    Function ret;
    ret.assignNode(External::create("./" + name + ".so", name));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  Function Function::external(const string& name, const string& bin_name,
                              const Dict& opts) {
    Function ret;
    ret.assignNode(External::create(bin_name, name));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  Function Function::external(const string& name, const Compiler& compiler,
                              const Dict& opts) {
    Function ret;
    ret.assignNode(External::create(compiler, name));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  bool Function::has_ivpsol(const string& name) {
    return Ivpsol::hasPlugin(name);
  }

  void Function::load_ivpsol(const string& name) {
    Ivpsol::loadPlugin(name);
  }

  string Function::doc_ivpsol(const string& name) {
    return Ivpsol::getPlugin(name).doc;
  }

  bool Function::has_qpsol(const string& name) {
    return Qpsol::hasPlugin(name);
  }

  void Function::load_qpsol(const string& name) {
    Qpsol::loadPlugin(name);
  }

  string Function::doc_qpsol(const string& name) {
    return Qpsol::getPlugin(name).doc;
  }

  bool Function::has_nlpsol(const string& name) {
    return Nlpsol::hasPlugin(name);
  }

  void Function::load_nlpsol(const string& name) {
    Nlpsol::loadPlugin(name);
  }

  string Function::doc_nlpsol(const string& name) {
    return Nlpsol::getPlugin(name).doc;
  }

  bool Function::has_nlsol(const string& name) {
    return Nlsol::hasPlugin(name);
  }

  void Function::load_nlsol(const string& name) {
    Nlsol::loadPlugin(name);
  }

  string Function::doc_nlsol(const string& name) {
    return Nlsol::getPlugin(name).doc;
  }

  bool Function::has_linsol(const string& name) {
    return Linsol::hasPlugin(name);
  }

  void Function::load_linsol(const string& name) {
    Linsol::loadPlugin(name);
  }

  string Function::doc_linsol(const string& name) {
    return Linsol::getPlugin(name).doc;
  }

  Function Function::nlsol_fun() {
    casadi_assert(!isNull());
    Nlsol* n = dynamic_cast<Nlsol*>(get());
    casadi_assert_message(n!=0, "Not a nlsol");
    return n->f_;
  }

  Function Function::nlsol_jac() {
    casadi_assert(!isNull());
    Nlsol* n = dynamic_cast<Nlsol*>(get());
    casadi_assert_message(n!=0, "Not a nlsol");
    return n->jac_;
  }

  Function Function::nlsol_linsol() {
    casadi_assert(!isNull());
    Nlsol* n = dynamic_cast<Nlsol*>(get());
    casadi_assert_message(n!=0, "Not a nlsol");
    return n->linsol_;
  }

  Function Function::ivpsol(const string& name, const string& solver,
                            const SXDict& dae, const Dict& opts) {
    return ivpsol(name, solver, Ivpsol::map2problem(dae), opts);
  }

  Function Function::ivpsol(const string& name, const string& solver,
                            const MXDict& dae, const Dict& opts) {
    return ivpsol(name, solver, Ivpsol::map2problem(dae), opts);
  }

  Function Function::ivpsol(const string& name, const string& solver,
                            const Function& dae, const Dict& opts) {
    if (dae.is_a("sxfunction")) {
      SXProblem p = Ivpsol::fun2problem<SX>(dae);
      return Function::ivpsol(name, solver, p, opts);
    } else {
      MXProblem p = Ivpsol::fun2problem<MX>(dae);
      return Function::ivpsol(name, solver, p, opts);
    }
  }

  Function Function::ivpsol(const string& name, const string& solver,
                                const pair<Function, Function>& dae,
                                const Dict& opts) {
    if (dae.first.is_a("sxfunction")) {
      SXProblem p = Ivpsol::fun2problem<SX>(dae.first, dae.second);
      return Function::ivpsol(name, solver, p, opts);
    } else {
      MXProblem p = Ivpsol::fun2problem<MX>(dae.first, dae.second);
      return Function::ivpsol(name, solver, p, opts);
    }
  }

  Function Function::ivpsol(const string& name, const string& solver,
                                const XProblem& dae, const Dict& opts) {
    Function ret;
    ret.assignNode(Ivpsol::getPlugin(solver).creator(name, dae));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  Function Function::ivpsol_dae() {
    casadi_assert(!isNull());
    Ivpsol* n = dynamic_cast<Ivpsol*>(get());
    casadi_assert_message(n!=0, "Not an integrator");
    return n->f_;
  }

  vector<string> Function::ivpsol_in() {
    vector<string> ret(ivpsol_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=ivpsol_in(i);
    return ret;
  }

  vector<string> Function::ivpsol_out() {
    vector<string> ret(ivpsol_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=ivpsol_out(i);
    return ret;
  }

  string Function::ivpsol_in(int ind) {
    switch (static_cast<IvpsolInput>(ind)) {
    case IVPSOL_X0:  return "x0";
    case IVPSOL_P:   return "p";
    case IVPSOL_Z0:  return "z0";
    case IVPSOL_RX0: return "rx0";
    case IVPSOL_RP:  return "rp";
    case IVPSOL_RZ0: return "rz0";
    case IVPSOL_NUM_IN: break;
    }
    return string();
  }

  string Function::ivpsol_out(int ind) {
    switch (static_cast<IvpsolOutput>(ind)) {
    case IVPSOL_XF:  return "xf";
    case IVPSOL_QF:  return "qf";
    case IVPSOL_ZF:  return "zf";
    case IVPSOL_RXF: return "rxf";
    case IVPSOL_RQF: return "rqf";
    case IVPSOL_RZF: return "rzf";
    case IVPSOL_NUM_OUT: break;
    }
    return string();
  }

  int Function::ivpsol_n_in() {
    return IVPSOL_NUM_IN;
  }

  int Function::ivpsol_n_out() {
    return IVPSOL_NUM_OUT;
  }

  Function Function::nlpsol(const string& name, const string& solver,
                                const SXDict& nlp, const Dict& opts) {
    return nlpsol(name, solver, Nlpsol::map2problem(nlp), opts);
  }

  Function Function::nlpsol(const string& name, const string& solver,
                                const MXDict& nlp, const Dict& opts) {
    return nlpsol(name, solver, Nlpsol::map2problem(nlp), opts);
  }

  Function Function::nlpsol(const string& name, const string& solver,
                                const Function& nlp, const Dict& opts) {
    if (nlp.is_a("sxfunction")) {
      return Function::nlpsol(name, solver,
                                  Nlpsol::fun2problem<SX>(nlp), opts);
    } else {
      return Function::nlpsol(name, solver,
                                  Nlpsol::fun2problem<MX>(nlp), opts);
    }
  }

  Function Function::nlpsol(const string& name, const string& solver,
                                const XProblem& nlp, const Dict& opts) {
    Function ret;
    ret.assignNode(Nlpsol::instantiatePlugin(name, solver, nlp));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  Function Function::nlpsol_nlp() {
    casadi_assert(!isNull());
    Nlpsol* n = dynamic_cast<Nlpsol*>(get());
    casadi_assert_message(n!=0, "Not an NLP solver");
    return n->nlp_;
  }

  Function Function::nlpsol_gradf() {
    casadi_assert(!isNull());
    Nlpsol* n = dynamic_cast<Nlpsol*>(get());
    casadi_assert_message(n!=0, "Not an NLP solver");
    return n->gradF();
  }

  Function Function::nlpsol_jacg() {
    casadi_assert(!isNull());
    Nlpsol* n = dynamic_cast<Nlpsol*>(get());
    casadi_assert_message(n!=0, "Not an NLP solver");
    return n->jacG();
  }

  Function Function::nlpsol_hesslag() {
    casadi_assert(!isNull());
    Nlpsol* n = dynamic_cast<Nlpsol*>(get());
    casadi_assert_message(n!=0, "Not an NLP solver");
    return n->hessLag();
  }

  vector<string> Function::nlpsol_in() {
    vector<string> ret(nlpsol_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=nlpsol_in(i);
    return ret;
  }

  vector<string> Function::nlpsol_out() {
    vector<string> ret(nlpsol_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=nlpsol_out(i);
    return ret;
  }

  string Function::nlpsol_in(int ind) {
    switch (static_cast<NlpsolInput>(ind)) {
    case NLPSOL_X0:     return "x0";
    case NLPSOL_P:      return "p";
    case NLPSOL_LBX:    return "lbx";
    case NLPSOL_UBX:    return "ubx";
    case NLPSOL_LBG:    return "lbg";
    case NLPSOL_UBG:    return "ubg";
    case NLPSOL_LAM_X0: return "lam_x0";
    case NLPSOL_LAM_G0: return "lam_g0";
    case NLPSOL_NUM_IN: break;
    }
    return string();
  }

  string Function::nlpsol_out(int ind) {
    switch (static_cast<NlpsolOutput>(ind)) {
    case NLPSOL_X:     return "x";
    case NLPSOL_F:     return "f";
    case NLPSOL_G:     return "g";
    case NLPSOL_LAM_X: return "lam_x";
    case NLPSOL_LAM_G: return "lam_g";
    case NLPSOL_LAM_P: return "lam_p";
    case NLPSOL_NUM_OUT: break;
    }
    return string();
  }

  int Function::nlpsol_n_in() {
    return NLPSOL_NUM_IN;
  }

  int Function::nlpsol_n_out() {
    return NLPSOL_NUM_OUT;
  }

  Function Function::qpsol(const string& name, const string& solver,
                               const SpDict& qp, const Dict& opts) {
    Function ret;
    ret.assignNode(Qpsol::instantiatePlugin(name, solver, qp));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  void Function::qpsol_debug(const string &filename) const {
    ofstream file;
    file.open(filename.c_str());
    qpsol_debug(file);
  }

  void Function::qpsol_debug(ostream &file) const {
    casadi_assert(!isNull());
    const Qpsol* n = dynamic_cast<const Qpsol*>(get());
    casadi_assert_message(n!=0, "Not a QP solver");
    return n->generateNativeCode(file);
  }

  vector<string> Function::qpsol_in() {
    vector<string> ret(qpsol_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=qpsol_in(i);
    return ret;
  }

  vector<string> Function::qpsol_out() {
    vector<string> ret(qpsol_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=qpsol_out(i);
    return ret;
  }

  string Function::qpsol_in(int ind) {
    switch (static_cast<QpsolInput>(ind)) {
    case QPSOL_H:      return "h";
    case QPSOL_G:      return "g";
    case QPSOL_A:      return "a";
    case QPSOL_LBA:    return "lba";
    case QPSOL_UBA:    return "uba";
    case QPSOL_LBX:    return "lbx";
    case QPSOL_UBX:    return "ubx";
    case QPSOL_X0:     return "x0";
    case QPSOL_LAM_X0: return "lam_x0";
    case QPSOL_NUM_IN: break;
    }
    return string();
  }

  string Function::qpsol_out(int ind) {
    switch (static_cast<QpsolOutput>(ind)) {
    case QPSOL_X:     return "x";
    case QPSOL_COST:  return "cost";
    case QPSOL_LAM_A: return "lam_a";
    case QPSOL_LAM_X: return "lam_x";
    case QPSOL_NUM_OUT: break;
    }
    return string();
  }

  int Function::qpsol_n_in() {
    return QPSOL_NUM_IN;
  }

  int Function::qpsol_n_out() {
    return QPSOL_NUM_OUT;
  }

  Function Function::linsol(const std::string& name, const std::string& solver,
                            const Sparsity& sp, int nrhs, const Dict& opts) {
    Function ret;
    if (solver=="none") {
      ret.assignNode(new Linsol(name, sp, nrhs));
    } else {
      ret.assignNode(Linsol::getPlugin(solver).creator(name, sp, nrhs));
    }
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  void Function::linsol_prepare(void* mem, const double** arg, double** res,
                                int* iw, double* w) {
    (*this)->linsol_prepare(mem, arg, res, iw, w);
  }

  void Function::linsol_solve(double* x, int nrhs, bool tr) {
    (*this)->linsol_solve(x, nrhs, tr);
  }

  void Function::linsol_solve(bool tr) {
    (*this)->linsol_solve(tr);
  }

  MX Function::linsol_solve(const MX& A, const MX& B, bool tr) {
    return (*this)->linsol_solve(A, B, tr);
  }

  void Function::linsol_solveL(double* x, int nrhs, bool tr) {
    (*this)->linsol_solveL(x, nrhs, tr);
  }

  Sparsity Function::linsol_cholesky_sparsity(bool tr) const {
    return (*this)->linsol_cholesky_sparsity(tr);
  }

  DMatrix Function::linsol_cholesky(bool tr) const {
    return (*this)->linsol_cholesky(tr);
  }

  void Function::linsol_spsolve(bvec_t* X, const bvec_t* B, bool tr) const {
    (*this)->linsol_spsolve(X, B, tr);
  }

  void Function::linsol_spsolve(DMatrix& X, const DMatrix& B, bool tr) const {
    (*this)->linsol_spsolve(X, B, tr);
  }

  Function Function::nlsol(const string& name, const string& solver,
                                const Dict& opts) const {
    Function ret;
    ret.assignNode(Nlsol::instantiatePlugin(name, solver, *this));
    ret.setOption(opts);
    ret.init();
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

  MemBlock::MemBlock() : mem_(0) {
  }

  MemBlock::MemBlock(const Function& f) : f_(f), mem_(0) {
    mem_ = f_.alloc_mem();
  }

  MemBlock::MemBlock(const MemBlock& obj) : f_(obj.f_), mem_(0) {
    mem_ = f_.alloc_mem();
  }

  MemBlock& MemBlock::operator=(const MemBlock& obj) {
    if (mem_) f_.free_mem(mem_);
    f_ = obj.f_;
    mem_ = f_.alloc_mem();
    return *this;
  }

  MemBlock::~MemBlock() {
    if (mem_) f_.free_mem(mem_);
  }

} // namespace casadi

