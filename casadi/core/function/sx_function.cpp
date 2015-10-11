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


#include "sx_function_internal.hpp"
#include "../std_vector_tools.hpp"
#include "../sx/sx_node.hpp"
#include "mx_function.hpp"

#include <limits>
#include <stack>
#include <deque>
#include <fstream>
#include <sstream>

namespace casadi {

  using namespace std;


  SXFunction::SXFunction() {
  }

  void SXFunction::construct(const std::string& name,
                             const std::vector<SX>& arg,
                             const std::vector<SX>& res,
                             const Dict& opts,
                             const std::vector<std::string>& ischeme,
                             const std::vector<std::string>& oscheme) {
    assignNode(new SXFunctionInternal(name, arg, res));
    if (!ischeme.empty()) setOption("input_scheme", ischeme);
    if (!oscheme.empty()) setOption("output_scheme", oscheme);
    setOption(opts);
    init();
  }

  SXFunction::SXFunction(const std::string& name, const Function &f, const Dict& opts) {
    vector<SX> arg = f.sx_in();
    vector<SX> res = Function(f)(arg);
    vector<string> name_in = f.name_in();
    vector<string> name_out = f.name_out();
    Dict opts2(opts);
    if (!name_in.empty() && !opts.count("input_scheme")) opts2["input_scheme"]=name_in;
    if (!name_out.empty() && !opts.count("output_scheme")) opts2["output_scheme"]=name_out;
    construct(name, arg, res, opts2);
  }

  SXFunction::SXFunction(const std::string& name, const std::vector<SX>& arg,
                         const std::vector<SX>& res, const Dict& opts) {
    construct(name, arg, res, opts);
  }

  SXFunction::SXFunction(const std::string& name, const pair<SXDict, vector<string> >& arg,
                         const std::vector<SX>& res, const Dict& opts) {
    construct(name, make_vector(arg), res, opts, arg.second);
  }

  SXFunction::SXFunction(const std::string& name, const std::vector<SX>& arg,
                         const pair<SXDict, vector<string> >& res, const Dict& opts) {
    construct(name, arg, make_vector(res), opts, std::vector<string>(), res.second);
  }

  SXFunction::SXFunction(const std::string& name, const pair<SXDict, vector<string> >& arg,
                         const pair<SXDict, vector<string> >& res, const Dict& opts) {
    construct(name, make_vector(arg), make_vector(res), opts, arg.second, res.second);
  }

#ifdef USE_CXX11
  SXFunction::SXFunction(const std::string& name, std::initializer_list<SX> arg,
                         std::initializer_list<SX> res, const Dict& opts) {
    construct(name, vector<SX>(arg), vector<SX>(res), opts);
  }

  SXFunction::SXFunction(const std::string& name, std::vector<SX> arg,
                         std::initializer_list<SX> res, const Dict& opts) {
    construct(name, arg, vector<SX>(res), opts);
  }

  SXFunction::SXFunction(const std::string& name, std::initializer_list<SX> arg,
                         std::vector<SX> res, const Dict& opts) {
    construct(name, vector<SX>(arg), res, opts);
  }
#endif // USE_CXX11

  const SXFunctionInternal* SXFunction::operator->() const {
    return static_cast<const SXFunctionInternal*>(Function::operator->());
  }

  SXFunctionInternal* SXFunction::operator->() {
    return static_cast<SXFunctionInternal*>(Function::operator->());
  }

  bool SXFunction::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const SXFunctionInternal*>(ptr)!=0;
  }

#ifdef WITH_DEPRECATED_FEATURES
  SX SXFunction::jac(int iind, int oind, bool compact, bool symmetric) {
    return (*this)->jac(iind, oind, compact, symmetric);
  }

  SX SXFunction::grad(int iind, int oind) {
    return (*this)->grad(iind, oind);
  }

  SX SXFunction::tang(int iind, int oind) {
    return (*this)->tang(iind, oind);
  }

  SX SXFunction::hess(int iind, int oind) {
    return (*this)->hess(iind, oind);
  }
#endif // WITH_DEPRECATED_FEATURES

  const vector<ScalarAtomic>& SXFunction::algorithm() const {
    return (*this)->algorithm_;
  }

  int SXFunction::countNodes() const {
    assertInit();
    return algorithm().size() - nnz_out();
  }

  void SXFunction::clearSymbolic() {
    (*this)->clearSymbolic();
  }

  SXFunction::SXFunction(const MXFunction& f) {
    SXFunction t("expand_" + f.name(), f);
    assignNode(t.get());
  }

  SXFunction::SXFunction(const Function& f) {
    const SXFunctionInternal* temp = dynamic_cast<const SXFunctionInternal*>(f.get());
    if (temp) {
      assignNode(const_cast<SXFunctionInternal*>(temp));
    } else {
      SXFunction t("expand_" + f.name(), f);
      assignNode(t.get());
    }
  }

  int SXFunction::getWorkSize() const {
    return (*this)->sz_w();
  }

} // namespace casadi
