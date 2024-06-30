/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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
#include "casadi_misc.hpp"
#include "sx_function.hpp"
#include "mx_function.hpp"
#include "switch.hpp"
#include "bspline.hpp"
#include "nlpsol.hpp"
#include "mapsum.hpp"
#include "conic.hpp"
#include "jit_function.hpp"
#include "serializing_stream.hpp"
#include "serializer.hpp"
#include "tools.hpp"

#include <cctype>
#include <fstream>
#include <typeinfo>

namespace casadi {
  // Throw informative error message
  #define THROW_ERROR(FNAME, WHAT) \
  throw CasadiException("Error in Function::" FNAME " for '" + this->name() + "' "\
    "[" + this->class_name() + "] at " + CASADI_WHERE + ":\n"\
    + std::string(WHAT));

  // Throw informative error message from constructor
  #define THROW_ERROR_NOOBJ(FNAME, WHAT, CLASS_NAME) \
  throw CasadiException("Error in Function::" FNAME " for '" + name + "' "\
      "[" CLASS_NAME "] at " + CASADI_WHERE + ":\n"\
      + std::string(WHAT));

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
    std::string tmp;
    while (true) {
      // Read a word
      std::streampos cur_pos = file.tellg();
      file >> tmp;
      if (!file.good()) return false;

      // Check if match
      if (str==tmp) return true;

      // If comment, continue to the end of the line
      if (tmp.at(0)=='#') {
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        continue;
      }

      // If mismatching name, rewind and break
      file.seekg(cur_pos);
      return false;
    }
  }

  Function::Function(const std::string& fname) {
    casadi_error("Not implemented");
  }

  Function::Function(const std::string& name,
                     const std::vector<SX>& ex_in, const std::vector<SX>& ex_out,
                     const Dict& opts) {
    construct(name, ex_in, ex_out, {}, {}, opts);
  }

  Function::Function(const std::string& name,
                     const std::vector<SX>& ex_in, const std::vector<SX>& ex_out,
                     const std::vector<std::string>& name_in,
                     const std::vector<std::string>& name_out,
                     const Dict& opts) {
    construct(name, ex_in, ex_out, name_in, name_out, opts);
  }

  Function::Function(const std::string& name,
                     const std::vector<MX>& ex_in, const std::vector<MX>& ex_out,
                     const Dict& opts) {
    construct(name, ex_in, ex_out, {}, {}, opts);
  }

  Function::Function(const std::string& name,
                     const std::vector<MX>& ex_in, const std::vector<MX>& ex_out,
                     const std::vector<std::string>& name_in,
                     const std::vector<std::string>& name_out,
                     const Dict& opts) {
    construct(name, ex_in, ex_out, name_in, name_out, opts);
  }

  Function::Function(const std::string& name,
      SXIList ex_in, const SXVector& ex_out, const Dict& opts) {
    construct(name, SXVector(ex_in), ex_out, {}, {}, opts);
  }

  Function::Function(const std::string& name,
      const SXVector& ex_in, SXIList ex_out, const Dict& opts) {
    construct(name, ex_in, SXVector(ex_out), {}, {}, opts);
  }

  Function::Function(const std::string& name, SXIList ex_in, SXIList ex_out, const Dict& opts) {
    construct(name, SXVector(ex_in), SXVector(ex_out), {}, {}, opts);
  }

  Function::Function(const std::string& name, SXIList ex_in, const SXVector& ex_out,
      const StringVector& name_in,
      const StringVector& name_out, const Dict& opts) {
    construct(name, SXVector(ex_in), ex_out, name_in, name_out, opts);
  }

  Function::Function(const std::string& name, const SXVector& ex_in, SXIList ex_out,
      const StringVector& name_in, const StringVector& name_out, const Dict& opts) {
    construct(name, ex_in, SXVector(ex_out), name_in, name_out, opts);
  }

  Function::Function(const std::string& name, SXIList ex_in, SXIList ex_out,
                     const StringVector& name_in, const StringVector& name_out, const Dict& opts) {
    construct(name, SXVector(ex_in), SXVector(ex_out), name_in, name_out, opts);
  }

  Function::Function(const std::string& name,
      MXIList ex_in, const MXVector& ex_out, const Dict& opts) {
    construct(name, MXVector(ex_in), ex_out, {}, {}, opts);
  }

  Function::Function(const std::string& name,
      const MXVector& ex_in, MXIList ex_out, const Dict& opts) {
    construct(name, ex_in, MXVector(ex_out), {}, {}, opts);
  }

  Function::Function(const std::string& name, MXIList ex_in, MXIList ex_out, const Dict& opts) {
    construct(name, MXVector(ex_in), MXVector(ex_out), {}, {}, opts);
  }

  Function::Function(const std::string& name, MXIList ex_in, const MXVector& ex_out,
      const StringVector& name_in, const StringVector& name_out, const Dict& opts) {
    construct(name, MXVector(ex_in), ex_out, name_in, name_out, opts);
  }

  Function::Function(const std::string& name, const MXVector& ex_in, MXIList ex_out,
      const StringVector& name_in, const StringVector& name_out, const Dict& opts) {
    construct(name, ex_in, MXVector(ex_out), name_in, name_out, opts);
  }

  Function::Function(const std::string& name, MXIList ex_in, MXIList ex_out,
      const StringVector& name_in, const StringVector& name_out, const Dict& opts) {
    construct(name, MXVector(ex_in), MXVector(ex_out), name_in, name_out, opts);
  }

  Function::Function(const std::string& name, const std::map<std::string, SX>& dict,
      const std::vector<std::string>& name_in, const std::vector<std::string>& name_out,
      const Dict& opts) {
    construct(name, dict, name_in, name_out, opts);
  }

  Function::Function(const std::string& name, const std::map<std::string, MX>& dict,
      const std::vector<std::string>& name_in, const std::vector<std::string>& name_out,
      const Dict& opts) {
    construct(name, dict, name_in, name_out, opts);
  }

  template<typename M>
  void Function::construct(const std::string& name, const std::map<std::string, M>& dict,
                           const std::vector<std::string>& name_in,
                           const std::vector<std::string>& name_out,
                           const Dict& opts) {
    std::vector<M> ex_in(name_in.size()), ex_out(name_out.size());
    for (auto&& i : dict) {
      std::vector<std::string>::const_iterator it;
      if ((it = std::find(name_in.begin(), name_in.end(), i.first))!=name_in.end()) {
        // Input expression
        ex_in[it-name_in.begin()] = i.second;
      } else if ((it = std::find(name_out.begin(), name_out.end(), i.first))!=name_out.end()) {
        // Output expression
        ex_out[it-name_out.begin()] = i.second;
      } else {
        // Neither
        casadi_error("Unknown dictionary entry: '" + i.first + "'");
      }
    }
    construct(name, ex_in, ex_out, name_in, name_out, opts);
  }

  void Function::construct(const std::string& name,
                           const std::vector<SX>& ex_in, const std::vector<SX>& ex_out,
                           const std::vector<std::string>& name_in,
                           const std::vector<std::string>& name_out,
                           const Dict& opts) {
    try {
      own(new SXFunction(name, ex_in, ex_out, name_in, name_out));
      (*this)->construct(opts);

      // Perform external transformations
      auto it = opts.find("external_transform");
      if (it!=opts.end()) {
        auto v = it->second.to_vector_vector();
        for (const std::vector<GenericType>& vec : v) {
          casadi_assert(vec.size()>=2, "external_transform: inner list must be length >=2");
          casadi_assert(vec.size()<=3, "external_transform: inner list must be length <=3");
          std::string name = vec[0].to_string();
          std::string op = vec[1].to_string();
          Dict opts = vec.size()==3 ? vec[2].to_dict() : Dict();
          operator=(external_transform(name, op, (*this), opts));
        }
      }

    } catch(std::exception& e) {
      THROW_ERROR_NOOBJ("Function", e.what(), "SXFunction");
    }
  }

  void Function::construct(const std::string& name,
                           const std::vector<MX>& ex_in, const std::vector<MX>& ex_out,
                           const std::vector<std::string>& name_in,
                           const std::vector<std::string>& name_out,
                           const Dict& opts) {
    try {
      own(new MXFunction(name, ex_in, ex_out, name_in, name_out));
      (*this)->construct(opts);

      // Perform post-construction expand
      auto it = opts.find("post_expand");
      if (it!=opts.end()) {
        if (!it->second) return;
        auto it = opts.find("post_expand_options");
        if (it==opts.end()) {
          operator=((*this).expand());
        } else {
          operator=((*this).expand((*this).name(), it->second));
        }
      }

      // Perform external transformations
      it = opts.find("external_transform");
      if (it!=opts.end()) {
        auto v = it->second.to_vector_vector();
        for (const std::vector<GenericType>& vec : v) {
          casadi_assert(vec.size()>=2, "external_transform: inner list must be length >=2");
          casadi_assert(vec.size()<=3, "external_transform: inner list must be length <=3");
          std::string name = vec[0].to_string();
          std::string op = vec[1].to_string();
          Dict opts = vec.size()==3 ? vec[2].to_dict() : Dict();
          operator=(external_transform(name, op, (*this), opts));
        }
      }

    } catch(std::exception& e) {
      THROW_ERROR_NOOBJ("Function", e.what(), "MXFunction");
    }
  }

  Function Function::jit(const std::string& name, const std::string& body,
                     const std::vector<std::string>& name_in,
                     const std::vector<std::string>& name_out,
                     const Dict& opts) {
    // Pass empty vectors -> default values
    std::vector<Sparsity> sparsity_in, sparsity_out;
    return jit(name, body, name_in, name_out, sparsity_in, sparsity_out, opts);
  }

  Function Function::jit(const std::string& name, const std::string& body,
                     const std::vector<std::string>& name_in,
                     const std::vector<std::string>& name_out,
                     const std::vector<Sparsity>& sparsity_in,
                     const std::vector<Sparsity>& sparsity_out,
                     const Dict& opts) {
    try {
      return create(new JitFunction(name, body, name_in, name_out,
                                    sparsity_in, sparsity_out), opts);
    } catch(std::exception& e) {
      THROW_ERROR_NOOBJ("jit", e.what(), "JitFunction");
    }
  }

  Function Function::expand() const {
    Dict opts;
    opts["ad_weight"] = (*this)->ad_weight();
    opts["ad_weight_sp"] = (*this)->sp_weight();
    opts["max_num_dir"] = (*this)->max_num_dir_;
    opts["is_diff_in"] = (*this)->is_diff_in_;
    opts["is_diff_out"] = (*this)->is_diff_out_;
    return expand(name(), opts);
  }

  Function Function::expand(const std::string& name, const Dict& opts) const {
    casadi_assert(!has_free(),
      "Function with free symbols cannot be expanded. "
      "List of free variables in your Function: " +
        join(get_free(), ","));

    Dict my_opts;
    my_opts["ad_weight"] = (*this)->ad_weight();
    my_opts["ad_weight_sp"] = (*this)->sp_weight();
    my_opts["max_num_dir"] = (*this)->max_num_dir_;
    my_opts["is_diff_in"] = (*this)->is_diff_in_;
    my_opts["is_diff_out"] = (*this)->is_diff_out_;
    update_dict(my_opts, opts);
    std::vector<SX> ex_in = sx_in();
    std::vector<SX> ex_out = Function(*this)(ex_in);
    return Function(name, ex_in, ex_out, name_in(), name_out(), my_opts);
  }

  Function Function::create(FunctionInternal* node) {
    Function ret;
    ret.own(node);
    return ret;
  }

  Function Function::create(FunctionInternal* node, const Dict& opts) {
    Function ret = create(node);
    ret->construct(opts);
    return ret;
  }

  FunctionInternal* Function::operator->() const {
    casadi_assert_dev(!is_null());
    return get();
  }

  FunctionInternal* Function::get() const {
    return static_cast<FunctionInternal*>(SharedObject::get());
  }

  void Function::call(const std::vector<DM> &arg, std::vector<DM> &res,
                      bool always_inline, bool never_inline) const {
    try {
      (*this)->call(arg, res, always_inline, never_inline);
    } catch(std::exception& e) {
      THROW_ERROR("call", e.what());
    }
  }

  void Function::call(const std::vector<SX> &arg, std::vector<SX>& res,
                      bool always_inline, bool never_inline) const {
    try {
      (*this)->call(arg, res, always_inline, never_inline);
    } catch(std::exception& e) {
      THROW_ERROR("call", e.what());
    }
  }

  void Function::call(const std::vector<MX> &arg, std::vector<MX>& res,
                      bool always_inline, bool never_inline) const {
    try {
      (*this)->call(arg, res, always_inline, never_inline);
    } catch(std::exception& e) {
      THROW_ERROR("call", e.what());
    }
  }

  std::vector<const double*> Function::buf_in(Function::VecArg arg) const {
    casadi_assert_dev(arg.size()==n_in());
    auto arg_it=arg.begin();
    std::vector<const double*> buf_arg(sz_arg());
    for (casadi_uint i=0; i<arg.size(); ++i) {
      casadi_assert_dev(arg_it->size()==nnz_in(i));
      buf_arg[i] = get_ptr(*arg_it++);
    }
    return buf_arg;
  }

  std::vector<double*> Function::buf_out(Function::VecRes res) const {
    res.resize(n_out());
    auto res_it=res.begin();
    std::vector<double*> buf_res(sz_res());
    for (casadi_uint i=0; i<res.size(); ++i) {
      res_it->resize(nnz_out(i));
      buf_res[i] = get_ptr(*res_it++);
    }
    return buf_res;
  }

  std::vector<double*> Function::buf_out(Function::VPrRes res) const {
    casadi_assert_dev(res.size()==n_out());
    auto res_it=res.begin();
    std::vector<double*> buf_res(sz_res());
    for (casadi_uint i=0; i<res.size(); ++i) {
      casadi_assert_dev(*res_it!=0);
      (*res_it)->resize(nnz_out(i));
      buf_res[i] = get_ptr(**res_it++);
    }
    return buf_res;
  }

  std::vector<const double*> Function::buf_in(Function::MapArg arg) const {
    // Return value (RVO)
    std::vector<const double*> ret(sz_arg(), nullptr);

    // Read inputs
    for (auto i=arg.begin(); i!=arg.end(); ++i) {
      casadi_int ind = index_in(i->first);
      casadi_assert_dev(i->second.size()==nnz_in(ind));
      ret[ind] = get_ptr(i->second);
    }

    return ret;
  }

  std::vector<double*> Function::buf_out(Function::MapRes res) const {
    // Return value (RVO)
    std::vector<double*> ret(sz_res(), nullptr);

    // Read outputs
    for (auto i=res.begin(); i!=res.end(); ++i) {
      casadi_int ind = index_out(i->first);
      i->second.resize(nnz_out(ind));
      ret[ind] = get_ptr(i->second);
    }

    return ret;
  }

  std::vector<double*> Function::buf_out(Function::MPrRes res) const {
    // Return value (RVO)
    std::vector<double*> ret(sz_res(), nullptr);

    // Read outputs
    for (auto i=res.begin(); i!=res.end(); ++i) {
      casadi_int ind = index_out(i->first);
      casadi_assert_dev(i->second!=0);
      i->second->resize(nnz_out(ind));
      ret[ind] = get_ptr(*i->second);
    }

    return ret;
  }

  template<typename D>
  void Function::call_gen(std::vector<const D*> arg, std::vector<D*> res) const {
    // Input buffer
    casadi_assert_dev(arg.size()>=n_in());
    arg.resize(sz_arg());

    // Output buffer
    casadi_assert_dev(res.size()>=n_out());
    res.resize(sz_res());

    // Work vectors
    std::vector<casadi_int> iw(sz_iw());
    std::vector<D> w(sz_w());

    // Evaluate memoryless
    (*this)(get_ptr(arg), get_ptr(res), get_ptr(iw), get_ptr(w), 0);
  }


  void Function::operator()(std::vector<const double*> arg, std::vector<double*> res) const {
    return call_gen(arg, res);
  }

  void Function::operator()(std::vector<const bvec_t*> arg, std::vector<bvec_t*> res) const {
    return call_gen(arg, res);
  }

  void Function::operator()(std::vector<const SXElem*> arg, std::vector<SXElem*> res) const {
    return call_gen(arg, res);
  }

  int Function::rev(std::vector<bvec_t*> arg, std::vector<bvec_t*> res) const {
    // Input buffer
    casadi_assert_dev(arg.size()>=n_in());
    arg.resize(sz_arg());

    // Output buffer
    casadi_assert_dev(res.size()>=n_out());
    res.resize(sz_res());

    // Work vectors
    std::vector<casadi_int> iw(sz_iw());
    std::vector<bvec_t> w(sz_w());

    // Evaluate memoryless
    return rev(get_ptr(arg), get_ptr(res), get_ptr(iw), get_ptr(w), 0);
  }

  Function Function::fold(casadi_int N, const Dict& opts) const {
    Function base = mapaccum(N, opts);
    std::vector<MX> base_in = base.mx_in();
    std::vector<MX> out = base(base_in);
    out[0] = out[0](Slice(), range((N-1)*size2_out(0), N*size2_out(0))); // NOLINT
    return Function("fold_"+name(), base_in, out, name_in(), name_out(), opts);
  }
  Function Function::mapaccum(casadi_int N, const Dict& opts) const {
    return mapaccum("mapaccum_"+name(), N, opts);
  }
  Function Function::mapaccum(const std::string& name, casadi_int N, const Dict& opts) const {
    return mapaccum(name, N, 1, opts);
  }
  Function Function::mapaccum(const std::string& name, casadi_int N, casadi_int n_accum,
                              const Dict& opts) const {
    Dict options = opts;

    // Default base
    casadi_int base = 10;
    auto it = options.find("base");
    if (it!=options.end()) {
      base = it->second;
      options.erase(it);
    }

    casadi_assert(N>0, "mapaccum: N must be positive");

    if (base==-1)
      return mapaccum(name, std::vector<Function>(N, *this), n_accum, options);
    casadi_assert(base>=2, "mapaccum: base must be positive");

    // Decompose N into
    std::vector<Function> chain;
    Function c = *this;
    while (N!=0) {
      casadi_int r = N % base;
      chain.insert(chain.end(), r, c);
      N = (N-r)/base;
      c = c.mapaccum(c.name()+"_acc"+str(base), std::vector<Function>(base, c), n_accum, options);
    }
    return mapaccum(name, chain, n_accum, options);
  }

  Function Function::mapaccum(const std::string& name,
                      const std::vector<Function>& chain, casadi_int n_accum,
                      const Dict& opts) const {
    // Shorthands
    casadi_int n_in = this->n_in(), n_out = this->n_out();
    // Consistency checks
    casadi_assert(!chain.empty(), "mapaccum: chain must be non-empty");
    casadi_assert(n_accum<=std::min(n_in, n_out), "mapaccum: too many accumulators");
    // Quick return?
    if (chain.size()==1) return chain[0];
    // Get symbolic expressions for inputs and outputs
    std::vector<MX> arg = mx_in();
    std::vector<MX> res;
    // Vectorized inputs and outputs
    std::vector<std::vector<MX>> varg(n_in), vres(n_out);
    for (casadi_int i=0; i<n_accum; ++i) varg[i].push_back(arg[i]);
    // For each function call
    for (const auto& f : chain) {

      // Stacked input expressions
      for (casadi_int i=n_accum; i<n_in; ++i) {
        arg[i] = MX::sym(name_in(i) + "_" + str(i), f.sparsity_in(i));
        varg[i].push_back(arg[i]);
      }

      // Call f
      res = f(arg);
      // Save output expressions
      for (casadi_int i=0; i<n_out; ++i) vres[i].push_back(res[i]);
      // Copy function output to input
      std::copy_n(res.begin(), n_accum, arg.begin());
      for (casadi_int i=0; i<n_accum; ++i) {
        // Ony get last component (allows nested calls)
        casadi_int ncol_out=f.size2_out(i), ncol_in=size2_in(i);
        if (ncol_out>ncol_in) {
          arg[i] = horzsplit(arg[i], {0, ncol_out-ncol_in, ncol_out}).back();
        }
      }
    }
    // Construct return
    for (casadi_int i=0; i<n_in; ++i) arg[i] = horzcat(varg[i]);
    for (casadi_int i=0; i<n_out; ++i) res[i] = horzcat(vres[i]);
    return Function(name, arg, res, name_in(), name_out(), opts);
  }

  Function Function::mapaccum(const std::string& name, casadi_int n,
                              const std::vector<casadi_int>& accum_in,
                              const std::vector<casadi_int>& accum_out,
                              const Dict& opts) const {
    // Shorthands
    casadi_int n_in = this->n_in(), n_out = this->n_out();
    // Consistency checks
    casadi_assert_dev(in_range(accum_in, n_in) && isUnique(accum_in));
    casadi_assert_dev(in_range(accum_out, n_out) && isUnique(accum_out));
    casadi_assert_dev(accum_in.size()==accum_out.size());
    casadi_int n_accum=accum_in.size();

    // Quick return if no need to reorder
    if (accum_in==range(n_accum) && accum_out==range(n_accum)) {
      return mapaccum(name, n, n_accum, opts);
    }

    // Need to do some reordering
    std::vector<casadi_int> temp_in = complement(accum_in, n_in);
    std::vector<casadi_int> order_in = accum_in;
    order_in.insert(order_in.end(), temp_in.begin(), temp_in.end());
    std::vector<casadi_int> temp_out = complement(accum_out, n_out);
    std::vector<casadi_int> order_out = accum_out;
    order_out.insert(order_out.end(), temp_out.begin(), temp_out.end());
    Function ret = slice("slice_" + name, order_in, order_out);
    ret = ret.mapaccum("mapacc_" + name, n, n_accum, opts);
    return ret.slice(name, lookupvector(order_in, n_in),
                     lookupvector(order_out, n_out), opts);
  }

  Function Function::mapaccum(const std::string& name, casadi_int n,
                              const std::vector<std::string>& accum_in,
                              const std::vector<std::string>& accum_out,
                              const Dict& opts) const {
    std::vector<casadi_int> accum_in_num, accum_out_num;
    for (const std::string& s : accum_in) accum_in_num.push_back(index_in(s));
    for (const std::string& s : accum_out) accum_out_num.push_back(index_out(s));
    return mapaccum(name, n, accum_in_num, accum_out_num, opts);
  }

  Function Function::map(casadi_int n,
    const std::vector<bool>& reduce_in,
    const std::vector<bool>& reduce_out,
    const Dict& opts) const {
    return MapSum::create("mapsum_" + str(n) + "_" + name(), "serial",
      *this, n, reduce_in, reduce_out, opts);
  }

  Function Function::map(const std::string& name, const std::string& parallelization, casadi_int n,
      const std::vector<casadi_int>& reduce_in, const std::vector<casadi_int>& reduce_out,
        const Dict& opts) const {
    // Wrap in an MXFunction
    Function f = map(n, parallelization);
    // Start with the fully mapped inputs
    std::vector<MX> arg = f.mx_in();
    std::vector<MX> f_arg = arg;
    // Replace reduced inputs
    for (casadi_int i : reduce_in) {
      arg[i] = mx_in(i);
      f_arg[i] = repmat(arg[i], 1, n);
    }
    // Get fully mapped outputs
    std::vector<MX> res = f(f_arg);
    // Replace reduced outputs
    for (casadi_int i : reduce_out) {
      res[i] = repsum(res[i], 1, n);
    }
    // Construct return
    return Function(name, arg, res, name_in(), name_out());
  }

  Function Function::map(const std::string& name, const std::string& parallelization, casadi_int n,
      const std::vector<std::string>& reduce_in, const std::vector<std::string>& reduce_out,
      const Dict& opts) const {
    std::vector<casadi_int> reduce_in_num, reduce_out_num;
    for (const std::string& s : reduce_in) reduce_in_num.push_back(index_in(s));
    for (const std::string& s : reduce_out) reduce_out_num.push_back(index_out(s));
    return map(name, parallelization, n, reduce_in_num, reduce_out_num, opts);
  }

  Function
  Function::map(casadi_int n, const std::string& parallelization,
      casadi_int max_num_threads) const {
    casadi_assert(max_num_threads>=1, "max_num_threads invalid.");
    // No need for logic when we are not saturating the limit
    if (n<=max_num_threads) return map(n, parallelization);

    // Floored division
    casadi_int d = n/max_num_threads;
    if (d*max_num_threads==n) {
      // Easy when n is divisable by max_num_threads
      return map(d, "serial").map(max_num_threads, parallelization);
    } else {
      // Create a base map that computes a bit too much
      Function base = map(d+1, "serial").map(max_num_threads, parallelization);
      std::vector<MX> ret_in, base_in;
      casadi_int rem = (d+1)*max_num_threads-n;
      for (casadi_int i=0;i<n_in();++i) {
        MX arg = MX::sym("arg", repmat(sparsity_in(i), 1, n));
        ret_in.push_back(arg);
        MX last_arg = arg(Slice(), range((n-1)*size2_in(i), n*size2_in(i))); // NOLINT
        base_in.push_back(horzcat(arg, repmat(last_arg, 1, rem)));
      }
      std::vector<MX> ret_out = base(base_in);
      for (casadi_int i=0;i<n_out();++i) {
        ret_out[i] = horzsplit(ret_out[i], {0, n*size2_out(i), ret_out[i].size2()})[0];
      }
      return Function("helper", ret_in, ret_out, name_in(), name_out());
    }
  }

  Function
  Function::map(casadi_int n, const std::string& parallelization) const {
    // Make sure not degenerate
    casadi_assert(n>0, "Degenerate map operation");
    // Quick return if possible
    if (n==1) return *this;
    // Unroll?
    if (parallelization=="unroll" || parallelization=="inline") {
      // Construct symbolic inputs
      std::vector<MX> arg(n_in());
      std::vector<std::vector<MX>> v(n, arg);
      std::vector<MX> tmp(n);
      for (casadi_int i=0; i<arg.size(); ++i) {
        for (casadi_int k=0; k<n; ++k) {
          tmp[k] = v[k][i] = MX::sym(name_in(i)+"_"+str(k), sparsity_in(i));
        }
        arg[i] = horzcat(tmp);
      }
      // Evaluate
      if (parallelization=="unroll") {
        for (auto&& w : v) w = (*this)(w);
      } else {
        for (auto&& w : v) call(std::vector<MX>(w), w, true, false);
      }
      // Gather outputs
      std::vector<MX> res(n_out());
      for (casadi_int i=0; i<res.size(); ++i) {
        for (casadi_int k=0; k<n; ++k) tmp[k] = v[k][i];
        res[i] = horzcat(tmp);
      }
      // Construct function
      return Function(name() + "_" + str(n), arg, res, name_in(), name_out());
    } else {
      // Generate/retrieve potentially cached map
      return (*this)->map(n, parallelization);
    }
  }

  Function Function::
  slice(const std::string& name, const std::vector<casadi_int>& order_in,
        const std::vector<casadi_int>& order_out, const Dict& opts) const {
    try {
      return (*this)->slice(name, order_in, order_out, opts);
    } catch(std::exception& e) {
      THROW_ERROR("slice", e.what());
    }
  }

  std::vector<MX> Function::mapsum(const std::vector< MX > &x,
                              const std::string& parallelization) const {
    try {
      return (*this)->mapsum_mx(x, parallelization);
    } catch(std::exception& e) {
      THROW_ERROR("mapsum", e.what());
    }
  }

  Function Function::conditional(const std::string& name, const std::vector<Function>& f,
                                 const Function& f_def, const Dict& opts) {
    try {
      return create(new Switch(name, f, f_def), opts);
    } catch(std::exception& e) {
      THROW_ERROR_NOOBJ("conditional", e.what(), "Switch");
    }
  }

  Function Function::conditional(const std::string& name,
      const Function& f, const Dict& opts) {
    try {
      // Create a dummy function with the same signature as f
      std::vector<MX> dummy_in = f.mx_in();
      std::vector<MX> dummy_out(f.n_out());
      for (casadi_int i = 0; i < dummy_out.size(); ++i) {
        dummy_out.at(i) = MX::zeros(f.sparsity_out(i));
      }
      Function dummy("dummy_" + f.name(), dummy_in, dummy_out, f.name_in(), f.name_out());
      // Form a conditional call
      return if_else(name, f, dummy, opts);
    } catch(std::exception& e) {
      THROW_ERROR_NOOBJ("conditional", e.what(), "Switch");
    }
  }

  Function Function::bspline(const std::string &name,
      const std::vector< std::vector<double> >& knots,
      const std::vector<double>& coeffs, const std::vector<casadi_int>& degree,
        casadi_int m, const Dict& opts) {
    try {
      casadi_assert(degree.size()==knots.size(), "Degree list length (" + str(degree.size()) + ") "
                  "must match knot list length (" + str(knots.size()) + ").");
      MX x = MX::sym("x", degree.size());
      std::vector<std::string> lookup_mode;
      Dict opts_remainder = extract_from_dict(opts, "lookup_mode", lookup_mode);
      Dict opts_bspline;
      opts_bspline["lookup_mode"] = lookup_mode;
      MX y = MX::bspline(x, DM(coeffs), knots, degree, m, opts_bspline);
      return Function(name, {x}, {y}, opts_remainder);
    } catch(std::exception& e) {
      THROW_ERROR_NOOBJ("bspline", e.what(), "BSpline");
    }
  }

  Function Function::if_else(const std::string& name, const Function& f_true,
                             const Function& f_false, const Dict& opts) {
    try {
      return create(new Switch(name, std::vector<Function>(1, f_false), f_true), opts);
    } catch(std::exception& e) {
      THROW_ERROR_NOOBJ("if_else", e.what(), "Switch");
    }
  }

  casadi_int Function::n_in() const {
    return (*this)->n_in_;
  }

  casadi_int Function::n_out() const {
    return (*this)->n_out_;
  }

  casadi_int Function::size1_in(casadi_int ind) const {
    return (*this)->size1_in(ind);
  }

  casadi_int Function::size2_in(casadi_int ind) const {
    return (*this)->size2_in(ind);
  }

  casadi_int Function::size1_out(casadi_int ind) const {
    return (*this)->size1_out(ind);
  }

  casadi_int Function::size2_out(casadi_int ind) const {
    return (*this)->size2_out(ind);
  }

  std::pair<casadi_int, casadi_int> Function::size_in(casadi_int ind) const {
    return (*this)->size_in(ind);
  }

  std::pair<casadi_int, casadi_int> Function::size_out(casadi_int ind) const {
    return (*this)->size_out(ind);
  }

  casadi_int Function::nnz_in() const {
    return (*this)->nnz_in();
  }

  casadi_int Function::nnz_out() const {
    return (*this)->nnz_out();
  }

  casadi_int Function::numel_in() const {
    return (*this)->numel_in();
  }

  casadi_int Function::numel_out() const {
    return (*this)->numel_out();
  }

  casadi_int Function::nnz_in(casadi_int ind) const {
    return (*this)->nnz_in(ind);
  }

  casadi_int Function::nnz_out(casadi_int ind) const {
    return (*this)->nnz_out(ind);
  }

  casadi_int Function::numel_in(casadi_int ind) const {
    return (*this)->numel_in(ind);
  }

  casadi_int Function::numel_out(casadi_int ind) const {
    return (*this)->numel_out(ind);
  }

  bool Function::uses_output() const {
    return (*this)->uses_output();
  }

#ifdef WITH_DEPRECATED_FEATURES
  Function Function::jacobian_old(casadi_int iind, casadi_int oind) const {
    // Redirect to factory class
    std::vector<std::string> s_in = name_in();
    std::vector<std::string> s_out = name_out();
    s_out.insert(s_out.begin(), "jac:" + name_out(oind) + ":" + name_in(iind));
    return factory(name() + "_jac", s_in, s_out);
  }

  Function Function::hessian_old(casadi_int iind, casadi_int oind) const {
    // Redirect to factory class
    std::vector<std::string> s_in = name_in();
    std::vector<std::string> s_out = name_out();
    s_out.insert(s_out.begin(), "grad:" + name_out(oind) + ":" + name_in(iind));
    s_out.insert(s_out.begin(),
                 "hess:" + name_out(oind) + ":" + name_in(iind) + ":" + name_in(iind));
    return factory(name() + "_hess", s_in, s_out);
  }

  const Sparsity Function::
  sparsity_jac(casadi_int iind, casadi_int oind, bool compact, bool symmetric) const {
    try {
      return (*this)->jac_sparsity(oind, iind, compact, symmetric);
    } catch(std::exception& e) {
      THROW_ERROR("sparsity_jac", e.what());
    }
  }
#endif  // WITH_DEPRECATED_FEATURES

  Function Function::jacobian() const {
    try {
      return (*this)->jacobian();
    } catch(std::exception& e) {
      THROW_ERROR("jacobian", e.what());
    }
  }

  bool Function::test_cast(const SharedObjectInternal* ptr) {
    return dynamic_cast<const FunctionInternal*>(ptr)!=nullptr;
  }

  Dict Function::stats(int mem) const {
    return (*this)->get_stats(memory(mem));
  }

  const std::vector<Sparsity>& Function::jac_sparsity(bool compact) const {
    // Make sure all are calculated
    for (casadi_int oind = 0; oind < n_out(); ++oind) {
      for (casadi_int iind = 0; iind < n_in(); ++iind) {
        (void)jac_sparsity(oind, iind, compact);
      }
    }
    // Return reference to internal cache
    return (*this)->jac_sparsity_[compact];
  }

  Sparsity Function::jac_sparsity(casadi_int oind, casadi_int iind, bool compact) const {
    try {
      return (*this)->jac_sparsity(oind, iind, compact, (*this)->jac_is_symm(oind, iind));
    } catch(std::exception& e) {
      THROW_ERROR("jac_sparsity", e.what());
    }
  }

  const std::vector<std::string>& Function::name_in() const {
    return (*this)->name_in_;
  }

  const std::vector<std::string>& Function::name_out() const {
    return (*this)->name_out_;
  }

  casadi_int Function::index_in(const std::string &name) const {
    try {
      return (*this)->index_in(name);
    } catch(std::exception& e) {
      THROW_ERROR("index_in", e.what());
    }
  }

  casadi_int Function::index_out(const std::string &name) const {
    try {
      return (*this)->index_out(name);
    } catch(std::exception& e) {
      THROW_ERROR("index_out", e.what());
    }
  }

  const std::string& Function::name_in(casadi_int ind) const {
    try {
      return (*this)->name_in_.at(ind);
    } catch(std::exception& e) {
      THROW_ERROR("name_in", e.what());
    }
  }

  const std::string& Function::name_out(casadi_int ind) const {
    try {
      return (*this)->name_out_.at(ind);
    } catch(std::exception& e) {
      THROW_ERROR("name_out", e.what());
    }
  }

  const Sparsity& Function::sparsity_in(casadi_int ind) const {
    try {
      return (*this)->sparsity_in_.at(ind);
    } catch(std::exception& e) {
      THROW_ERROR("sparsity_in", e.what());
    }
  }

  const Sparsity& Function::sparsity_in(const std::string &iname) const {
    try {
      return sparsity_in(index_in(iname));
    } catch(std::exception& e) {
      THROW_ERROR("sparsity_in", e.what());
    }
  }

  const Sparsity& Function::sparsity_out(casadi_int ind) const {
    try {
      return (*this)->sparsity_out_.at(ind);
    } catch(std::exception& e) {
      THROW_ERROR("sparsity_out", e.what());
    }
  }

  const Sparsity& Function::sparsity_out(const std::string &iname) const {
    try {
      return sparsity_out(index_out(iname));
    } catch(std::exception& e) {
      THROW_ERROR("sparsity_out", e.what());
    }
  }

  bool Function::is_diff_in(casadi_int ind) const {
    try {
      return (*this)->is_diff_in_.at(ind);
    } catch(std::exception& e) {
      THROW_ERROR("is_diff_in", e.what());
    }
  }

  bool Function::is_diff_out(casadi_int ind) const {
    try {
      return (*this)->is_diff_out_.at(ind);
    } catch(std::exception& e) {
      THROW_ERROR("is_diff_out", e.what());
    }
  }

  std::vector<bool> Function::is_diff_in() const {
    try {
      return (*this)->is_diff_in_;
    } catch(std::exception& e) {
      THROW_ERROR("is_diff_in", e.what());
    }
  }

  std::vector<bool> Function::is_diff_out() const {
    try {
      return (*this)->is_diff_out_;
    } catch(std::exception& e) {
      THROW_ERROR("is_diff_out", e.what());
    }
  }

  void Function::sz_work(size_t& sz_arg, size_t& sz_res, size_t& sz_iw, size_t& sz_w) const {
    (*this)->sz_work(sz_arg, sz_res, sz_iw, sz_w);
  }

  size_t Function::sz_arg() const { return (*this)->sz_arg();}

  size_t Function::sz_res() const { return (*this)->sz_res();}

  size_t Function::sz_iw() const { return (*this)->sz_iw();}

  size_t Function::sz_w() const { return (*this)->sz_w();}

  int Function::operator()(const bvec_t** arg, bvec_t** res,
                            casadi_int* iw, bvec_t* w, int mem) const {
    try {
      return (*this)->sp_forward(arg, res, iw, w, memory(mem));
    } catch(std::exception& e) {
      THROW_ERROR("operator()", e.what());
    }
  }

  int Function::rev(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, int mem) const {
    try {
      return (*this)->sp_reverse(arg, res, iw, w, memory(mem));
    } catch(std::exception& e) {
      THROW_ERROR("rev", e.what());
    }
  }

  void Function::set_work(const double**& arg, double**& res, casadi_int*& iw, double*& w,
                          int mem) const {
    try {
      (*this)->set_work(memory(mem), arg, res, iw, w);
    } catch(std::exception& e) {
      THROW_ERROR("set_work", e.what());
    }
  }

  void Function::set_temp(const double** arg, double** res, casadi_int* iw, double* w,
                          int mem) const {
    try {
      (*this)->set_temp(memory(mem), arg, res, iw, w);
    } catch(std::exception& e) {
      THROW_ERROR("set_temp", e.what());
    }
  }

  void Function::setup(const double** arg, double** res, casadi_int* iw, double* w,
                          int mem) const {
    try {
      (*this)->setup(memory(mem), arg, res, iw, w);
    } catch(std::exception& e) {
      THROW_ERROR("setup", e.what());
    }
  }

  Function Function::forward(casadi_int nfwd) const {
    try {
      return (*this)->forward(nfwd);
    } catch(std::exception& e) {
      THROW_ERROR("forward", e.what());
    }
  }

  Function Function::reverse(casadi_int nadj) const {
    try {
      return (*this)->reverse(nadj);
    } catch(std::exception& e) {
      THROW_ERROR("reverse", e.what());
    }
  }

  void Function::print_dimensions(std::ostream &stream) const {
    (*this)->print_dimensions(stream);
  }

  void Function::print_options(std::ostream &stream) const {
    (*this)->print_options(stream);
  }

  void Function::print_option(const std::string &name, std::ostream &stream) const {
    (*this)->print_option(name, stream);
  }

  bool Function::has_option(const std::string &option_name) const {
    try {
      return (*this)->has_option(option_name);
    } catch(std::exception& e) {
      THROW_ERROR("has_option", e.what());
      return false;  // never reached
    }
  }

  void Function::change_option(const std::string& option_name, const GenericType& option_value) {
    try {
      // Assert existance
      if (!has_option(option_name))
        casadi_error("Option '" + option_name + "' does not exist");
      // Call internal class
      (*this)->change_option(option_name, option_value);
    } catch(std::exception& e) {
      THROW_ERROR("change_option", e.what());
    }
  }

  std::vector<std::string> Function::get_free() const {
    return (*this)->get_free();
  }

  std::string Function::generate(const Dict& opts) const {
    return generate(name(), opts);
  }

  std::string Function::generate(const std::string& fname, const Dict& opts) const {
    return generate(fname, "", opts);
  }

  std::string Function::generate(const std::string& fname, const std::string& path, const Dict& opts) const {
    CodeGenerator gen(fname, opts);
    gen.add(*this);
    return gen.generate(path);
  }

  std::string Function::generate_dependencies(const std::string& fname, const Dict& opts) const {
    return (*this)->generate_dependencies(fname, opts);
  }

  void Function::generate_in(const std::string& fname, const std::vector<DM>& arg) {
    std::vector<double> d = nz_from_in(arg);

    // Set up output stream
    std::ofstream of(fname);
    casadi_assert(of.good(), "Error opening stream '" + fname + "'.");
    normalized_setup(of);

    // Encode each output
    for (casadi_int i=0; i<d.size(); ++i) {
      normalized_out(of, d[i]);
      of << std::endl;
    }
  }

  void Function::generate_out(const std::string& fname, const std::vector<DM>& res) {
    std::vector<double> d = nz_from_out(res);

    // Set up output stream
    std::ofstream of(fname);
    casadi_assert(of.good(), "Error opening stream '" + fname + "'.");
    normalized_setup(of);

    // Encode each output
    for (casadi_int i=0; i<d.size(); ++i) {
      normalized_out(of, d[i]);
      of << std::endl;
    }
  }

  std::vector<DM> Function::generate_in(const std::string& fname) {
    DM data = DM::from_file(fname, "txt");
    // Empty files are okay
    if (data.is_empty(true)) data = DM(0, 1);
    casadi_assert(data.is_vector() && data.is_dense(), "Expected dense vector");
    casadi_assert(data.numel()==nnz_in(),
      "Dimension mismatch: file contains a vector of size " + str(data.numel())
      + ", while size " + str(nnz_in()) + " was expected.");

    return nz_to_in(data.nonzeros());
  }

  std::vector<DM> Function::generate_out(const std::string& fname) {
    DM data = DM::from_file(fname, "txt");
    // Empty files are okay
    if (data.is_empty(true)) data = DM(0, 1);
    casadi_assert(data.is_vector() && data.is_dense(), "Expected dense vector");
    casadi_assert(data.numel()==nnz_out(),
      "Dimension mismatch: file contains a vector of size " + str(data.numel())
      + ", while size " + str(nnz_out()) + " was expected.");

    return nz_to_out(data.nonzeros());
  }

  void Function::export_code(const std::string& lang,
      std::ostream &stream, const Dict& options) const {
    return (*this)->export_code(lang, stream, options);
  }

  void Function::export_code(const std::string& lang,
      const std::string &fname, const Dict& options) const {
    std::ofstream stream(fname);
    return (*this)->export_code(lang, stream, options);
  }


  void Function::save(const std::string &fname, const Dict& opts) const {
    FileSerializer fs(fname, opts);
    fs.pack(*this);
  }

  std::string Function::serialize(const Dict& opts) const {
    std::stringstream ss;
    serialize(ss, opts);
    return ss.str();
  }

  void Function::serialize(std::ostream &stream, const Dict& opts) const {
    SerializingStream s(stream, opts);
    return serialize(s);
  }

  void Function::serialize(SerializingStream &s) const {
    if (is_null()) {
      s.pack("Function::null", true);
    } else {
      s.pack("Function::null", false);
      (*this)->serialize(s);
    }
  }

  Function Function::deserialize(DeserializingStream& s) {
    bool is_null;
    s.unpack("Function::null", is_null);
    if (is_null) return Function();
    return FunctionInternal::deserialize(s);
  }

  std::string Function::export_code(const std::string& lang, const Dict& options) const {
    std::stringstream ss;
    (*this)->export_code(lang, ss, options);
    return ss.str();
  }

  const std::string& Function::name() const {
    if (is_null()) {
      static std::string null = "null";
      return null;
    } else {
      return (*this)->name_;
    }
  }

  bool Function::check_name(const std::string& name) {
    // Check if empty
    if (name.empty()) return false;

    // Check if keyword
    for (const char* kw : {"null", "jac", "hess"}) {
      if (name==kw) return false;
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

  Function Function::deserialize(std::istream& stream) {
    DeserializingStream s(stream);
    return deserialize(s);
  }

  Function Function::load(const std::string& filename) {
    FileDeserializer fs(filename);
    auto t = fs.pop_type();
    if (t==SerializerBase::SerializationType::SERIALIZED_FUNCTION) {
      return fs.blind_unpack_function();
    } else {
      casadi_error("File is not loadable with 'load'. Use 'FileDeserializer' instead.");
    }
  }

  Function Function::deserialize(const std::string& s) {
    std::stringstream ss;
    ss << s;
    return deserialize(ss);
  }

  std::string Function::fix_name(const std::string& name) {
    // Quick return if already valid name
    if (check_name(name)) return name;

    // If empty, name it "unnamed"
    if (name.empty()) return "unnamed";

    // Construct a sane name
    std::stringstream ss;

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
      if (ss.str()==kw) ss << "1";
    }

    return ss.str();
  }

  std::vector<DM> Function::operator()(const std::vector<DM>& arg) const {
    std::vector<DM> res;
    call(arg, res);
    return res;
  }

  std::vector<SX> Function::operator()(const std::vector<SX>& arg) const {
    std::vector<SX> res;
    call(arg, res);
    return res;
  }

  std::vector<MX> Function::operator()(const std::vector<MX>& arg) const {
    std::vector<MX> res;
    call(arg, res);
    return res;
  }

  template<typename M>
  void Function::call_gen(const std::map<std::string, M>& arg, std::map<std::string, M>& res,
                       bool always_inline, bool never_inline) const {
    // Convert to vector arguments
    std::vector<M> arg_v = (*this)->convert_arg(arg);

    // Make call
    std::vector<M> res_v;
    call(arg_v, res_v, always_inline, never_inline);

    // Save to map
    res.clear();
    for (casadi_int i=0; i<res_v.size(); ++i) {
      res[name_out(i)] = res_v[i];
    }
  }

  const DMDict Function::operator()(const DMDict& arg) const {
    DMDict res;
    call(arg, res);
    return res;
  }

  const SXDict Function::operator()(const SXDict& arg) const {
    SXDict res;
    call(arg, res);
    return res;
  }

  const MXDict Function::operator()(const MXDict& arg) const {
    MXDict res;
    call(arg, res);
    return res;
  }

  void Function::call(const DMDict& arg, DMDict& res,
                      bool always_inline, bool never_inline) const {
    try {
      call_gen(arg, res, always_inline, never_inline);
    } catch(std::exception& e) {
      THROW_ERROR("call", e.what());
    }
  }

  void Function::call(const SXDict& arg, SXDict& res,
                      bool always_inline, bool never_inline) const {
    try {
      call_gen(arg, res, always_inline, never_inline);
    } catch(std::exception& e) {
      THROW_ERROR("call", e.what());
    }
  }

  void Function::call(const MXDict& arg, MXDict& res,
                      bool always_inline, bool never_inline) const {
    try {
      call_gen(arg, res, always_inline, never_inline);
    } catch(std::exception& e) {
      THROW_ERROR("call", e.what());
    }
  }

  double Function::default_in(casadi_int ind) const {
    return (*this)->get_default_in(ind);
  }

  double Function::max_in(casadi_int ind) const {
    return (*this)->get_max_in(ind);
  }

  double Function::min_in(casadi_int ind) const {
    return (*this)->get_min_in(ind);
  }

  std::vector<double> Function::nominal_in(casadi_int ind) const {
    return (*this)->get_nominal_in(ind);
  }

  std::vector<double> Function::nominal_out(casadi_int ind) const {
    return (*this)->get_nominal_out(ind);
  }

#ifdef WITH_EXTRA_CHECKS
  // Initialize at zero depth
  thread_local casadi_int Function::call_depth_ = 0;
#endif // WITH_EXTRA_CHECKS

  int Function::operator()(const double** arg, double** res,
      casadi_int* iw, double* w) const {
    scoped_checkout<Function> mem(*this);
    return operator()(arg, res, iw, w, mem);
  }

  int Function::operator()(const double** arg, double** res,
      casadi_int* iw, double* w, int mem) const {
    try {
#ifdef WITH_EXTRA_CHECKS
      // Should never happen
      casadi_assert_dev(call_depth_>=0);
      call_depth_++;
      // For consistency check
      casadi_int depth = call_depth_;
#endif // WITH_EXTRA_CHECKS
      int ret = (*this)->eval_gen(arg, res, iw, w, memory(mem));
#ifdef WITH_EXTRA_CHECKS
      // Consitency check
      casadi_assert_dev(call_depth_==depth);
      call_depth_--;
#endif // WITH_EXTRA_CHECKS
      return ret;
    } catch (KeyboardInterruptException& e) {
      (void)e;  // unused
#ifdef WITH_EXTRA_CHECKS
      call_depth_--;
#endif // WITH_EXTRA_CHECKS
      throw;
    } catch(std::exception& e) {
#ifdef WITH_EXTRA_CHECKS
      call_depth_--;
#endif // WITH_EXTRA_CHECKS
      (*this)->print_in(uerr(), arg, true);
      THROW_ERROR("operator()", e.what());
    }
  }

  int Function::operator()(const SXElem** arg, SXElem** res,
      casadi_int* iw, SXElem* w, int mem) const {
    try {
      return (*this)->eval_sx(arg, res, iw, w, memory(mem));
    } catch(std::exception& e) {
      THROW_ERROR("operator()", e.what());
    }
  }

  const SX Function::sx_in(casadi_int iind) const {
    try {
      return (*this)->sx_in(iind);
    } catch(std::exception& e) {
      THROW_ERROR("sx_in", e.what());
    }
  }

  const SX Function::sx_out(casadi_int oind) const {
    try {
      return (*this)->sx_out(oind);
    } catch(std::exception& e) {
      THROW_ERROR("sx_out", e.what());
    }
  }

  const std::vector<SX> Function::sx_in() const {
    try {
      return (*this)->sx_in();
    } catch(std::exception& e) {
      THROW_ERROR("sx_in", e.what());
    }
  }

  const std::vector<SX> Function::sx_out() const {
    try {
      return (*this)->sx_out();
    } catch(std::exception& e) {
      THROW_ERROR("sx_out", e.what());
    }
  }

  const MX Function::mx_in(casadi_int ind) const {
    return (*this)->mx_in(ind);
  }

  const MX Function::mx_out(casadi_int ind) const {
    return (*this)->mx_out(ind);
  }

  const std::vector<MX> Function::mx_in() const {
    return (*this)->mx_in();
  }

  const std::vector<MX> Function::mx_out() const {
    return (*this)->mx_out();
  }

  std::vector<double> Function::nz_from_in(const std::vector<DM>& arg) const {
    return (*this)->nz_in(arg);
  }

  std::vector<double> Function::nz_from_out(const std::vector<DM>& res) const {
    return (*this)->nz_out(res);
  }

  std::vector<DM> Function::nz_to_in(const std::vector<double>& arg) const {
    return (*this)->nz_in(arg);
  }

  std::vector<DM> Function::nz_to_out(const std::vector<double>& res) const {
    return (*this)->nz_out(res);
  }

  DMDict Function::convert_in(const std::vector<DM>& arg) const {
    return (*this)->convert_arg(arg);
  }

  std::vector<DM> Function::convert_in(const DMDict& arg) const {
    return (*this)->convert_arg(arg);
  }

  DMDict Function::convert_out(const std::vector<DM>& arg) const {
    return (*this)->convert_res(arg);
  }

  std::vector<DM> Function::convert_out(const DMDict& arg) const {
    return (*this)->convert_res(arg);
  }

  SXDict Function::convert_in(const std::vector<SX>& arg) const {
    return (*this)->convert_arg(arg);
  }

  std::vector<SX> Function::convert_in(const SXDict& arg) const {
    return (*this)->convert_arg(arg);
  }

  SXDict Function::convert_out(const std::vector<SX>& arg) const {
    return (*this)->convert_res(arg);
  }

  std::vector<SX> Function::convert_out(const SXDict& arg) const {
    return (*this)->convert_res(arg);
  }

  MXDict Function::convert_in(const std::vector<MX>& arg) const {
    return (*this)->convert_arg(arg);
  }

  std::vector<MX> Function::convert_in(const MXDict& arg) const {
    return (*this)->convert_arg(arg);
  }

  MXDict Function::convert_out(const std::vector<MX>& arg) const {
    return (*this)->convert_res(arg);
  }

  std::vector<MX> Function::convert_out(const MXDict& arg) const {
    return (*this)->convert_res(arg);
  }

  bool Function::is_a(const std::string& type, bool recursive) const {
    return (*this)->is_a(type, recursive);
  }

  std::vector<SX> Function::free_sx() const {
    try {
      return (*this)->free_sx();
    } catch(std::exception& e) {
      THROW_ERROR("free_sx", e.what());
    }
  }

  std::vector<MX> Function::free_mx() const {
    try {
      return (*this)->free_mx();
    } catch(std::exception& e) {
      THROW_ERROR("free_mx", e.what());
    }
  }

  bool Function::has_spfwd() const {
    return (*this)->has_spfwd();
  }

  bool Function::has_sprev() const {
    return (*this)->has_sprev();
  }

  bool Function::has_free() const {
    return (*this)->has_free();
  }

  void Function::generate_lifted(Function& vdef_fcn, Function& vinit_fcn) const {
    try {
      (*this)->generate_lifted(vdef_fcn, vinit_fcn);
    } catch(std::exception& e) {
      THROW_ERROR("generate_lifted", e.what());
    }
  }

  casadi_int Function::n_instructions() const {
    try {
      return (*this)->n_instructions();
    } catch(std::exception& e) {
      THROW_ERROR("n_instructions", e.what());
    }
  }

  MX Function::instruction_MX(casadi_int k) const {
    try {
      return (*this)->instruction_MX(k);
    } catch(std::exception& e) {
      THROW_ERROR("instruction_MX", e.what());
    }
  }

  SX Function::instructions_sx() const {
    try {
      return (*this)->instructions_sx();
    } catch(std::exception& e) {
      THROW_ERROR("instructions_sx", e.what());
    }
  }

  casadi_int Function::instruction_id(casadi_int k) const {
    try {
      return (*this)->instruction_id(k);
    } catch(std::exception& e) {
      THROW_ERROR("instruction_id", e.what());
    }
  }

  std::vector<casadi_int> Function::instruction_input(casadi_int k) const {
    try {
      return (*this)->instruction_input(k);
    } catch(std::exception& e) {
      THROW_ERROR("instruction_input", e.what());
    }
  }

  double Function::instruction_constant(casadi_int k) const {
    try {
      return (*this)->instruction_constant(k);
    } catch(std::exception& e) {
      THROW_ERROR("instruction_constant", e.what());
    }
  }

  std::vector<casadi_int> Function::instruction_output(casadi_int k) const {
    try {
      return (*this)->instruction_output(k);
    } catch(std::exception& e) {
      THROW_ERROR("instruction_output", e.what());
    }
  }

  casadi_int Function::n_nodes() const {
    try {
      return (*this)->n_nodes();
    } catch(std::exception& e) {
      THROW_ERROR("n_nodes", e.what());
    }
  }

  casadi_int Function::checkout() const {
    return (*this)->checkout();
  }

  void Function::release(int mem) const {
    (*this)->release(mem);
  }

  void* Function::memory(int ind) const {
    return (*this)->memory(ind);
  }

  void Function::assert_size_in(casadi_int i, casadi_int nrow, casadi_int ncol) const {
    casadi_assert(size1_in(i)==nrow && size2_in(i)==ncol,
      "Incorrect shape for " + str(*this) + " input " + str(i) + " \""
      + name_in(i) + "\". Expected " + str(nrow) + "-by-" + str(ncol)
      + " but got " + str(size1_in(i)) +  "-by-" + str(size2_in(i)));

  }

  void Function::assert_size_out(casadi_int i, casadi_int nrow, casadi_int ncol) const {
    casadi_assert(size1_out(i)==nrow && size2_out(i)==ncol,
      "Incorrect shape for " + str(*this) + " output " + str(i) + " \""
      + name_out(i) + "\". Expected " + str(nrow) + "-by-" + str(ncol)
      + " but got " + str(size1_out(i)) +  "-by-" + str(size2_out(i)));
  }

  void Function::assert_sparsity_out(casadi_int i, const Sparsity& sp,
      casadi_int n, bool allow_all_zero_sparse) const {
    // Assert shape
    assert_size_out(i, sp.size1(), sp.size2() * n);
    // Quick return if empty sparse
    if (allow_all_zero_sparse && sparsity_out(i).nnz() == 0) return;
    // Check sparsities
    casadi_assert(sparsity_out(i).is_stacked(sp, n), "Mismatching sparsity "
      "(but correct dimensions) for " + str(*this) + " output " + name_out(i));
  }

  Function Function::
  factory(const std::string& name,
          const std::vector<std::string>& s_in,
          const std::vector<std::string>& s_out,
          const AuxOut& aux,
          const Dict& opts) const {
     try {
       return (*this)->factory(name, s_in, s_out, aux, opts);
     } catch(std::exception& e) {
       THROW_ERROR("factory", "Failed to create " + name + ":" + str(s_in) + "->" + str(s_out)
        + " with " + str(aux) + ":\n" + str(e.what()));
     }
  }

  std::vector<bool> Function::
  which_depends(const std::string& s_in, const std::vector<std::string>& s_out,
      casadi_int order, bool tr) const {
    try {
      return (*this)->which_depends(s_in, s_out, order, tr);
    } catch(std::exception& e) {
      THROW_ERROR("which_depends", e.what());
    }
  }

  Dict Function::cache() const {
    try {
      return (*this)->cache();
    } catch(std::exception& e) {
      THROW_ERROR("cache", e.what());
      return {};
    }
  }

  std::vector<std::string> Function::get_function() const {
    try {
      return (*this)->get_function();
    } catch(std::exception& e) {
      THROW_ERROR("get_function", e.what());
      return {};
    }
  }

  Function Function::get_function(const std::string &name) const {
    try {
      return (*this)->get_function(name);
    } catch(std::exception& e) {
      THROW_ERROR("get_function", e.what());
    }
  }

  bool Function::has_function(const std::string& fname) const {
    try {
      return (*this)->has_function(fname);
    } catch(std::exception& e) {
      THROW_ERROR("has_function", e.what());
      return false;
    }
  }

  std::vector<Function> Function::find_functions(casadi_int max_depth) const {
    try {
      // If negative, make largest positive number
      if (max_depth < 0) max_depth = std::numeric_limits<casadi_int>::max();
      // The internal routine uses a map to avoid duplicate functions
      std::map<FunctionInternal*, Function> all_fun;
      (*this)->find(all_fun, max_depth);
      // Create return object
      std::vector<Function> ret;
      ret.reserve(all_fun.size());
      for (auto&& e : all_fun) ret.push_back(e.second);
      return ret;
    } catch(std::exception& e) {
      THROW_ERROR("find", e.what());
      return {};
    }
  }

  Function Function::find_function(const std::string &name, casadi_int max_depth) const {
    try {
      // If negative, make largest positive number
      if (max_depth < 0) max_depth = std::numeric_limits<casadi_int>::max();
      // The internal routine uses a map to avoid duplicate functions
      std::map<FunctionInternal*, Function> all_fun;
      (*this)->find(all_fun, max_depth);
      // Search this map
      for (auto&& e : all_fun) {
        if (e.second.name() == name) return e.second;
      }
      // Not found
      casadi_error("'" + name + "' not found");
    } catch(std::exception& e) {
      THROW_ERROR("find", e.what());
      return Function();
    }
  }


  Function Function::oracle() const {
    try {
      return (*this)->oracle();
    } catch(std::exception& e) {
      THROW_ERROR("oracle", e.what());
    }
  }

  Function Function::wrap() const {
    return (*this)->wrap();
  }

  Function Function::wrap_as_needed(const Dict& opts) const {
    return (*this)->wrap_as_needed(opts);
  }

  bool Function::operator==(const Function& f) const {
    try {
      casadi_assert(!is_null(), "lhs is null");
      casadi_assert(!f.is_null(), "rhs is null");
      return get()==f.get();
    } catch(std::exception& e) {
      THROW_ERROR("operator==", e.what());
    }
  }

  Dict Function::info() const {
    return (*this)->info();
  }

  FunctionBuffer::FunctionBuffer(const Function& f) : f_(f) {
    w_.resize(f_.sz_w());
    iw_.resize(f_.sz_iw());
    arg_.resize(f_.sz_arg());
    res_.resize(f_.sz_res());
    mem_ = f_->checkout();
    mem_internal_ = f.memory(mem_);
    f_node_ = f.operator->();
  }

  FunctionBuffer::~FunctionBuffer() {
    if (f_->release_) {
      f_->release_(mem_);
    } else {
      f_.release(mem_);
    }
  }

  FunctionBuffer::FunctionBuffer(const FunctionBuffer& f) : f_(f.f_) {
    operator=(f);
  }

  FunctionBuffer& FunctionBuffer::operator=(const FunctionBuffer& f) {
    f_ = f.f_;
    w_ = f.w_; iw_ = f.iw_; arg_ = f.arg_; res_ = f.res_; f_node_ = f.f_node_;
    // Checkout fresh memory
    if (f_->checkout_) {
      mem_ = f_->checkout_();
    } else {
      mem_ = f_.checkout();
      mem_internal_ = f_.memory(mem_);
    }

    return *this;
  }

  void FunctionBuffer::set_arg(casadi_int i, const double* a, casadi_int size) {
    casadi_assert(size>=f_.nnz_in(i)*sizeof(double),
     "Buffer is not large enough. Needed " + str(f_.nnz_in(i)*sizeof(double)) +
     " bytes, got " + str(size) + ".");
    arg_.at(i) = a;
  }
  void FunctionBuffer::set_res(casadi_int i, double* a, casadi_int size) {
    casadi_assert(size>=f_.nnz_out(i)*sizeof(double),
     "Buffer is not large enough. Needed " + str(f_.nnz_out(i)*sizeof(double)) +
     " bytes, got " + str(size) + ".");
    res_.at(i) = a;
  }
  void FunctionBuffer::_eval() {
    if (f_node_->eval_) {
      ret_ = f_node_->eval_(get_ptr(arg_), get_ptr(res_), get_ptr(iw_), get_ptr(w_), mem_);
    } else {
      ret_ = f_node_->eval(get_ptr(arg_), get_ptr(res_), get_ptr(iw_), get_ptr(w_), mem_internal_);
    }
  }
  int FunctionBuffer::ret() {
    return ret_;
  }

  void CASADI_EXPORT _function_buffer_eval(void* raw) {
    static_cast<FunctionBuffer*>(raw)->_eval();
  }

} // namespace casadi
