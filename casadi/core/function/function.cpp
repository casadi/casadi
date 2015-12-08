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
#include "integrator.hpp"
#include "qpsol.hpp"
#include "nlpsol.hpp"
#include "rootfinder.hpp"
#include "linsol.hpp"
#include "jit.hpp"

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
    vector<const double*> ret(sz_arg(), 0);

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
    vector<const double*> ret(sz_arg(), 0);

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
    (*this)(getPtr(arg), getPtr(res), getPtr(buf_iw), getPtr(buf_w), 0);
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
    // Get pointers to input arguments
    int n_in = this->n_in();
    vector<const double*> arg(sz_arg());
    for (int i=0; i<n_in; ++i) arg[i]=input(i).ptr();

    // Get pointers to output arguments
    int n_out = this->n_out();
    vector<double*> res(sz_res());
    for (int i=0; i<n_out; ++i) res[i]=output(i).ptr();

    // Temporary memory
    std::vector<int> iw(sz_iw());
    std::vector<double> w(sz_w());

    // Call memory-less
    (*this)->eval(getPtr(arg), getPtr(res), getPtr(iw), getPtr(w), 0);
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

  size_t Function::sz_iw() const { return (*this)->sz_iw();}

  size_t Function::sz_w() const { return (*this)->sz_w();}

  void* Function::alloc_mem() { return (*this)->alloc_mem();}

  void Function::free_mem(void* mem) { (*this)->free_mem(mem);}

  void Function::operator()(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, void* mem) {
    (*this)->spFwd(arg, res, iw, w, mem);
  }

  void Function::rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, void* mem) {
    (*this)->spAdj(arg, res, iw, w, mem);
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
    if (isNull()) {
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

  void Function::operator()(const double** arg, double** res, int* iw, double* w, void* mem) {
    (*this)->eval(arg, res, iw, w, mem);
  }

  void Function::operator()(const SXElem** arg, SXElem** res, int* iw, SXElem* w, void* mem) {
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

  void Function::init() {
    (*this)->init();
    (*this)->finalize();
  }

  Function external(const string& name, const Dict& opts) {
    Function ret;
    ret.assignNode(External::create("./" + name + ".so", name));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  Function external(const string& name, const string& bin_name,
                              const Dict& opts) {
    Function ret;
    ret.assignNode(External::create(bin_name, name));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  Function external(const string& name, const Compiler& compiler,
                              const Dict& opts) {
    Function ret;
    ret.assignNode(External::create(compiler, name));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  bool has_integrator(const string& name) {
    return Integrator::hasPlugin(name);
  }

  void load_integrator(const string& name) {
    Integrator::loadPlugin(name);
  }

  string doc_integrator(const string& name) {
    return Integrator::getPlugin(name).doc;
  }

  bool has_qpsol(const string& name) {
    return Qpsol::hasPlugin(name);
  }

  void load_qpsol(const string& name) {
    Qpsol::loadPlugin(name);
  }

  string doc_qpsol(const string& name) {
    return Qpsol::getPlugin(name).doc;
  }

  bool has_nlpsol(const string& name) {
    return Nlpsol::hasPlugin(name);
  }

  void load_nlpsol(const string& name) {
    Nlpsol::loadPlugin(name);
  }

  string doc_nlpsol(const string& name) {
    return Nlpsol::getPlugin(name).doc;
  }

  bool has_rootfinder(const string& name) {
    return Rootfinder::hasPlugin(name);
  }

  void load_rootfinder(const string& name) {
    Rootfinder::loadPlugin(name);
  }

  string doc_rootfinder(const string& name) {
    return Rootfinder::getPlugin(name).doc;
  }

  bool has_linsol(const string& name) {
    return Linsol::hasPlugin(name);
  }

  void load_linsol(const string& name) {
    Linsol::loadPlugin(name);
  }

  string doc_linsol(const string& name) {
    return Linsol::getPlugin(name).doc;
  }

  Function Function::rootfinder_fun() {
    casadi_assert(!isNull());
    Rootfinder* n = dynamic_cast<Rootfinder*>(get());
    casadi_assert_message(n!=0, "Not a rootfinder");
    return n->f_;
  }

  Function Function::rootfinder_jac() {
    casadi_assert(!isNull());
    Rootfinder* n = dynamic_cast<Rootfinder*>(get());
    casadi_assert_message(n!=0, "Not a rootfinder");
    return n->jac_;
  }

  Function Function::rootfinder_linsol() {
    casadi_assert(!isNull());
    Rootfinder* n = dynamic_cast<Rootfinder*>(get());
    casadi_assert_message(n!=0, "Not a rootfinder");
    return n->linsol_;
  }

  Function integrator(const string& name, const string& solver,
                            const SXDict& dae, const Dict& opts) {
    return integrator(name, solver, Integrator::map2problem(dae), opts);
  }

  Function integrator(const string& name, const string& solver,
                            const MXDict& dae, const Dict& opts) {
    return integrator(name, solver, Integrator::map2problem(dae), opts);
  }

  Function integrator(const string& name, const string& solver,
                            const Function& dae, const Dict& opts) {
    if (dae.is_a("sxfunction")) {
      SXProblem p = Integrator::fun2problem<SX>(dae);
      return integrator(name, solver, p, opts);
    } else {
      MXProblem p = Integrator::fun2problem<MX>(dae);
      return integrator(name, solver, p, opts);
    }
  }

  Function integrator(const string& name, const string& solver,
                                const pair<Function, Function>& dae,
                                const Dict& opts) {
    if (dae.first.is_a("sxfunction")) {
      SXProblem p = Integrator::fun2problem<SX>(dae.first, dae.second);
      return integrator(name, solver, p, opts);
    } else {
      MXProblem p = Integrator::fun2problem<MX>(dae.first, dae.second);
      return integrator(name, solver, p, opts);
    }
  }

  Function integrator(const string& name, const string& solver,
                                const XProblem& dae, const Dict& opts) {
    Function ret;
    ret.assignNode(Integrator::getPlugin(solver).creator(name, dae));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

  Function Function::integrator_dae() {
    casadi_assert(!isNull());
    Integrator* n = dynamic_cast<Integrator*>(get());
    casadi_assert_message(n!=0, "Not an integrator");
    return n->f_;
  }

  vector<string> integrator_in() {
    vector<string> ret(integrator_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=integrator_in(i);
    return ret;
  }

  vector<string> integrator_out() {
    vector<string> ret(integrator_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=integrator_out(i);
    return ret;
  }

  string integrator_in(int ind) {
    switch (static_cast<IntegratorInput>(ind)) {
    case INTEGRATOR_X0:  return "x0";
    case INTEGRATOR_P:   return "p";
    case INTEGRATOR_Z0:  return "z0";
    case INTEGRATOR_RX0: return "rx0";
    case INTEGRATOR_RP:  return "rp";
    case INTEGRATOR_RZ0: return "rz0";
    case INTEGRATOR_NUM_IN: break;
    }
    return string();
  }

  string integrator_out(int ind) {
    switch (static_cast<IntegratorOutput>(ind)) {
    case INTEGRATOR_XF:  return "xf";
    case INTEGRATOR_QF:  return "qf";
    case INTEGRATOR_ZF:  return "zf";
    case INTEGRATOR_RXF: return "rxf";
    case INTEGRATOR_RQF: return "rqf";
    case INTEGRATOR_RZF: return "rzf";
    case INTEGRATOR_NUM_OUT: break;
    }
    return string();
  }

  int integrator_n_in() {
    return INTEGRATOR_NUM_IN;
  }

  int integrator_n_out() {
    return INTEGRATOR_NUM_OUT;
  }

  Function nlpsol(const string& name, const string& solver,
                                const SXDict& nlp, const Dict& opts) {
    return nlpsol(name, solver, Nlpsol::map2problem(nlp), opts);
  }

  Function nlpsol(const string& name, const string& solver,
                                const MXDict& nlp, const Dict& opts) {
    return nlpsol(name, solver, Nlpsol::map2problem(nlp), opts);
  }

  Function nlpsol(const string& name, const string& solver,
                                const Function& nlp, const Dict& opts) {
    if (nlp.is_a("sxfunction")) {
      return nlpsol(name, solver, Nlpsol::fun2problem<SX>(nlp), opts);
    } else {
      return nlpsol(name, solver, Nlpsol::fun2problem<MX>(nlp), opts);
    }
  }

  Function nlpsol(const string& name, const string& solver,
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

  vector<string> nlpsol_in() {
    vector<string> ret(nlpsol_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=nlpsol_in(i);
    return ret;
  }

  vector<string> nlpsol_out() {
    vector<string> ret(nlpsol_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=nlpsol_out(i);
    return ret;
  }

  string nlpsol_in(int ind) {
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

  string nlpsol_out(int ind) {
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

  int nlpsol_n_in() {
    return NLPSOL_NUM_IN;
  }

  int nlpsol_n_out() {
    return NLPSOL_NUM_OUT;
  }

  Function qpsol(const string& name, const string& solver,
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

  vector<string> qpsol_in() {
    vector<string> ret(qpsol_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=qpsol_in(i);
    return ret;
  }

  vector<string> qpsol_out() {
    vector<string> ret(qpsol_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=qpsol_out(i);
    return ret;
  }

  string qpsol_in(int ind) {
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

  string qpsol_out(int ind) {
    switch (static_cast<QpsolOutput>(ind)) {
    case QPSOL_X:     return "x";
    case QPSOL_COST:  return "cost";
    case QPSOL_LAM_A: return "lam_a";
    case QPSOL_LAM_X: return "lam_x";
    case QPSOL_NUM_OUT: break;
    }
    return string();
  }

  int qpsol_n_in() {
    return QPSOL_NUM_IN;
  }

  int qpsol_n_out() {
    return QPSOL_NUM_OUT;
  }

  Function linsol(const std::string& name, const std::string& solver,
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

  bool Function::has_mem() const {
    return (*this)->has_mem();
  }

  void Function::linsol_factorize(Memory& m, const double* A) {
    (*this)->linsol_factorize(m, A);
  }

  void Function::linsol_solve(Memory& m, double* x, int nrhs, bool tr) {
    (*this)->linsol_solve(m, x, nrhs, tr);
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

  DM Function::linsol_cholesky(bool tr) const {
    return (*this)->linsol_cholesky(tr);
  }

  void Function::linsol_spsolve(bvec_t* X, const bvec_t* B, bool tr) const {
    (*this)->linsol_spsolve(X, B, tr);
  }

  void Function::linsol_spsolve(DM& X, const DM& B, bool tr) const {
    (*this)->linsol_spsolve(X, B, tr);
  }

  Function rootfinder(const std::string& name, const std::string& solver,
                   const Function& f, const Dict& opts) {
    Function ret;
    ret.assignNode(Rootfinder::instantiatePlugin(name, solver, f));
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

  Memory::Memory()
    : f(0), arg(0), res(0), iw(0), w(0), mem(0), own_(false) {
  }

  Memory::Memory(FunctionInternal *_f, const double** _arg, double** _res,
                 int* _iw, double* _w, void* _mem)
    : f(_f), arg(_arg), res(_res), iw(_iw), w(_w), mem(_mem), own_(false) {
    setup(arg, res, iw, w);
  }

  Memory::Memory(const Function& _f, const double** _arg, double** _res,
                 int* _iw, double* _w, void* _mem)
    : f(_f.get()), arg(_arg), res(_res), iw(_iw), w(_w), mem(_mem), own_(false) {
    setup(arg, res, iw, w);
  }

  Memory::Memory(const Function& _f)
    : f(_f.get()), arg(0), res(0), iw(0), w(0), mem(0), own_(true) {
    // Make an owning reference
    casadi_assert(f!=0);
    f->count++; // prevent object from being deleted

    // Allocate work vectors
    arg = new const double*[f->sz_arg()];
    res = new double*[f->sz_res()];
    iw = new int[f->sz_iw()];
    w = new double[f->sz_w()];

    // Allocate memory
    mem = f->alloc_mem();

    // Set up memory object
    setup(arg, res, iw, w);
  }

  Memory::Memory(Memory&& obj)
    : f(obj.f), arg(obj.arg), res(obj.res), iw(obj.iw), w(obj.w), mem(obj.mem),
        own_(obj.own_) {
  }

  Memory::~Memory() {
    if (own_) {
      // Free work vectors
      if (arg) delete[] arg;
      if (res) delete[] res;
      if (iw) delete[] iw;
      if (w) delete[] w;

      // Free memory object
      if (mem) f->free_mem(mem);
      if (--f->count==0) delete f;
      f = 0;
    }
  }

  void Memory::setup(const double** arg1, double** res1, int* iw1, double* w1) {
    casadi_assert(f!=0);
    arg1 += f->n_in();
    res1 += f->n_out();
    f->set_work(*this, arg1, res1, iw1, w1);
    f->set_temp(*this, arg1, res1, iw1, w1);
  }

  Memory& Memory::operator=(Memory&& obj) {
    if (&obj != this) {
      // Copy the members
      f = obj.f;
      arg = obj.arg;
      res = obj.res;
      iw = obj.iw;
      w = obj.w;
      mem = obj.mem;
      own_ = obj.own_;
      // Clear from assigning object
      obj.f = 0;
      obj.arg = 0;
      obj.res = 0;
      obj.iw = 0;
      obj.w = 0;
      obj.mem = 0;
      obj.own_ = 0;
    }
    return *this;
  }

  template<typename M>
  Function qpsol_nlp(const std::string& name, const std::string& solver,
                     const Problem<M>& qp, const Dict& opts) {
    // We have: minimize    f(x) = 1/2 * x' H x + c'x
    //          subject to  lbx <= x <= ubx
    //                      lbg <= g(x) = A x + b <= ubg
    M x = qp.in[NL_X];
    M p = qp.in[NL_P];
    M f = qp.out[NL_F];
    M g = qp.out[NL_G];
    if (g.is_empty(true)) g = M(0, 1); // workaround

    // Gradient of the objective: gf == Hx + g
    M gf = M::gradient(f, x);

    // Identify the linear term in the objective
    M c = substitute(gf, x, M::zeros(x.sparsity()));

    // Identify the quadratic term in the objective
    M H = M::jacobian(gf, x, true);

    // Identify the constant term in the constraints
    M b = substitute(g, x, M::zeros(x.sparsity()));

    // Identify the linear term in the constraints
    M A = M::jacobian(g, x);

    // Create a function for calculating the required matrices vectors
    Function prob(name + "_qp", qp.in, {H, c, A, b});

    // Create the QP solver
    Function qpsol_f = qpsol(name + "_qpsol", solver,
                             {{"h", H.sparsity()}, {"a", A.sparsity()}}, opts);

    // Create an MXFunction with the right signature
    vector<MX> ret_in(NLPSOL_NUM_IN);
    ret_in[NLPSOL_X0] = MX::sym("x0", x.sparsity());
    ret_in[NLPSOL_P] = MX::sym("p", p.sparsity());
    ret_in[NLPSOL_LBX] = MX::sym("lbx", x.sparsity());
    ret_in[NLPSOL_UBX] = MX::sym("ubx", x.sparsity());
    ret_in[NLPSOL_LBG] = MX::sym("lbg", g.sparsity());
    ret_in[NLPSOL_UBG] = MX::sym("ubg", g.sparsity());
    ret_in[NLPSOL_LAM_X0] = MX::sym("lam_x0", x.sparsity());
    ret_in[NLPSOL_LAM_G0] = MX::sym("lam_g0", g.sparsity());
    vector<MX> ret_out(NLPSOL_NUM_OUT);

    // Get expressions for the QP matrices and vectors
    vector<MX> v(NL_NUM_IN);
    v[NL_X] = ret_in[NLPSOL_X0];
    v[NL_P] = ret_in[NLPSOL_P];
    v = prob(v);

    // Call the QP solver
    vector<MX> w(QPSOL_NUM_IN);
    w[QPSOL_H] = v.at(0);
    w[QPSOL_G] = v.at(1);
    w[QPSOL_A] = v.at(2);
    w[QPSOL_LBX] = ret_in[NLPSOL_LBX];
    w[QPSOL_UBX] = ret_in[NLPSOL_UBX];
    w[QPSOL_LBA] = ret_in[NLPSOL_LBG] - v.at(3);
    w[QPSOL_UBA] = ret_in[NLPSOL_UBG] - v.at(3);
    w[QPSOL_X0] = ret_in[NLPSOL_X0];
    w[QPSOL_LAM_X0] = ret_in[NLPSOL_LAM_X0];
    w = qpsol_f(w);

    // Get expressions for the solution
    ret_out[NLPSOL_X] = w[QPSOL_X];
    ret_out[NLPSOL_F] = w[QPSOL_COST];
    ret_out[NLPSOL_G] = mul(v.at(2), w[QPSOL_X]) + v.at(3);
    ret_out[NLPSOL_LAM_X] = w[QPSOL_LAM_X];
    ret_out[NLPSOL_LAM_G] = w[QPSOL_LAM_A];
    ret_out[NLPSOL_LAM_P] = MX::nan(p.sparsity());
    return Function(name, ret_in, ret_out, nlpsol_in(), nlpsol_out());
  }

  Function qpsol(const std::string& name, const std::string& solver,
                 const XProblem& qp, const Dict& opts) {
    if (qp.is_sx) {
      return qpsol_nlp<SX>(name, solver, qp, opts);
    } else {
      return qpsol_nlp<MX>(name, solver, qp, opts);
    }
  }

  Function qpsol(const std::string& name, const std::string& solver,
                 const SXDict& qp, const Dict& opts) {
    return qpsol(name, solver, Nlpsol::map2problem(qp), opts);
  }

  Function qpsol(const std::string& name, const std::string& solver,
                 const MXDict& qp, const Dict& opts) {
    return qpsol(name, solver, Nlpsol::map2problem(qp), opts);
  }

  Function jit(const std::string& name, int n_in, int n_out,
               const std::string& body, const Dict& opts) {
    Function ret;
    ret.assignNode(new Jit(name, n_in, n_out, body, opts));
    ret.setOption(opts);
    ret.init();
    return ret;
  }

} // namespace casadi

