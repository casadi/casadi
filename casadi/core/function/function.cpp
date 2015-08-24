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
#include "../function/mx_function.hpp"
#include <typeinfo>
#include "../std_vector_tools.hpp"
#include "../function/map.hpp"
#include "../function/mapaccum.hpp"


using namespace std;

namespace casadi {

  Function::Function() {
  }

  Function::~Function() {
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

  Function Function::mapaccum(const std::string& name, int N, const Dict& options) const {
    std::vector<bool> accum_input(nIn(), false);
    accum_input[0] = true;
    std::vector<int> accum_output(1, 0);
    return MapAccum(name, *this, N, accum_input, accum_output, false, options);
  }

  Function Function::map(const std::string& name, int N, const Dict& options) const {
    std::vector<bool> repeated_input(nIn(), true);
    std::vector<bool> repeated_output(nOut(), true);
    return Map(name, *this, N, repeated_input, repeated_output, options);
  }

  vector<vector<MX> > Function::map(const vector<vector<MX> > &x,
                                    const std::string& parallelization) {
    assertInit();
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
                                    const std::string& parallelization) {
    assertInit();
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

    std::vector<bool> repeat_n;
    for (int i=0;i<x.size();++i) {
      repeat_n.push_back(x[i].size2()/input(i).size2()==n);
    }

    // Call the internal function

    Dict options;
    options["parallelization"] = parallelization;
    Function ms = Map("map", *this, n, repeat_n, std::vector<bool>(nOut(), true), options);

    // Call the internal function
    return ms(x);
  }

  vector<MX> Function::mapsum(const vector< MX > &x,
                                    const std::string& parallelization) {
    assertInit();
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

    std::vector<bool> repeat_n;
    for (int i=0;i<x.size();++i) {
      repeat_n.push_back(x[i].size2()/input(i).size2()==n);
    }

    Dict options;
    options["parallelization"] = parallelization;
    Function ms = Map("mapsum", *this, n, repeat_n, std::vector<bool>(nOut(), false), options);

    // Call the internal function
    return ms(x);
  }

  void Function::evaluate() {
    assertInit();
    (*this)->evaluate();
  }

  int Function::nIn() const {
    return (*this)->nIn();
  }

  int Function::nOut() const {
    return (*this)->nOut();
  }

  int Function::nnzIn() const {
    return (*this)->nnzIn();
  }

  int Function::nnzOut() const {
    return (*this)->nnzOut();
  }

  int Function::numelIn() const {
    return (*this)->numelIn();
  }

  int Function::numelOut() const {
    return (*this)->numelOut();
  }

  Function Function::jacobian(int iind, int oind, bool compact, bool symmetric) {
    assertInit();
    return (*this)->jacobian(iind, oind, compact, symmetric);
  }

  void Function::setJacobian(const Function& jac, int iind, int oind, bool compact) {
    (*this)->setJacobian(jac, iind, oind, compact);
  }

  Function Function::gradient(int iind, int oind) {
    assertInit();
    return (*this)->gradient(iind, oind);
  }

  Function Function::tangent(int iind, int oind) {
    assertInit();
    return (*this)->tangent(iind, oind);
  }

  Function Function::hessian(int iind, int oind) {
    assertInit();
    return (*this)->hessian(iind, oind);
  }

  Function Function::fullJacobian() {
    assertInit();
    return (*this)->fullJacobian();
  }

  void Function::setFullJacobian(const Function& jac) {
    (*this)->full_jacobian_ = jac;
  }

  bool Function::testCast(const SharedObjectNode* ptr) {
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

  const Sparsity Function::jacSparsity(int iind, int oind, bool compact, bool symmetric) {
    return (*this)->jacSparsity(iind, oind, compact, symmetric);
  }

  void Function::setJacSparsity(const Sparsity& sp, int iind, int oind, bool compact) {
    (*this)->setJacSparsity(sp, iind, oind, compact);
  }

  std::vector<MX> Function::symbolicInput(bool unique) const {
    if (unique) {
      return (*this)->FunctionInternal::symbolicInput();
    } else {
      return (*this)->symbolicInput();
    }
  }

  std::vector<SX> Function::symbolicInputSX() const {
    return (*this)->symbolicInputSX();
  }

  std::vector<MX> Function::symbolicOutput() const {
    return (*this)->symbolicOutput();
  }

  std::vector<std::string> Function::inputScheme() const {
    return (*this)->ischeme_;
  }

  std::vector<std::string> Function::outputScheme() const {
    return (*this)->oscheme_;
  }

  int Function::inputIndex(const std::string &name) const {
    return (*this)->inputIndex(name);
  }

  int Function::outputIndex(const std::string &name) const {
    return (*this)->outputIndex(name);
  }

  std::string Function::inputName(int ind) const {
    return (*this)->inputName(ind);
  }

  std::string Function::outputName(int ind) const {
    return (*this)->outputName(ind);
  }

  std::string Function::inputDescription(int ind) const {
    return (*this)->inputDescription(ind);
  }

  std::string Function::outputDescription(int ind) const {
    return (*this)->outputDescription(ind);
  }

  const Matrix<double>& Function::input(int i) const {
    return (*this)->input(i);
  }

  const Matrix<double>& Function::input(const std::string &iname) const {
    return (*this)->input(iname);
  }

  Sparsity Function::inputSparsity(int ind) const {
    return (*this)->inputSparsity(ind);
  }

  Sparsity Function::inputSparsity(const std::string &iname) const {
    return (*this)->inputSparsity(iname);
  }

  Sparsity Function::outputSparsity(int ind) const {
    return (*this)->outputSparsity(ind);
  }

  Sparsity Function::outputSparsity(const std::string &iname) const {
    return (*this)->outputSparsity(iname);
  }

  Matrix<double>& Function::input(int i) {
    return (*this)->input(i);
  }

  Matrix<double>& Function::input(const std::string &iname) {
    return (*this)->input(iname);
  }

  const Matrix<double>& Function::output(int i) const {
    return (*this)->output(i);
  }

  const Matrix<double>& Function::output(const std::string &oname) const {
    return (*this)->output(oname);
  }

  Matrix<double>& Function::output(int i) {
    return (*this)->output(i);
  }

  Matrix<double>& Function::output(const std::string &oname) {
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

  void Function::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    (*this)->spFwdSwitch(arg, res, iw, w);
  }

  void Function::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) {
    (*this)->spAdjSwitch(arg, res, iw, w);
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
    vector<MX> arg = symbolicInput();
    vector<MX> res = (*this)(arg);
    vector<MX> ret_in(arg), ret_out(res);

    // Number inputs and outputs
    int num_in = nIn();
    int num_out = nOut();

    // Forward sensitivities
    if (nfwd>0) {
      Function dfcn = derForward(nfwd);
      arg = dfcn.symbolicInput();
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
      Function dfcn = derReverse(nadj);
      arg = dfcn.symbolicInput();
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
    ss << "derivative_" << getOption("name") << "_" << nfwd << "_" << nadj;

    // Names of inputs
    std::vector<std::string> i_names;
    i_names.reserve(nIn()*(1+nfwd)+nOut()*nadj);
    const std::vector<std::string>& ischeme=(*this)->ischeme_;
    const std::vector<std::string>& oscheme=(*this)->oscheme_;

    // Nondifferentiated inputs
    for (int i=0; i<nIn(); ++i) {
      i_names.push_back("der_" + ischeme.at(i));
    }

    // Forward seeds
    for (int d=0; d<nfwd; ++d) {
      for (int i=0; i<nIn(); ++i) {
        ss.str(string());
        ss << "fwd" << d << "_" << ischeme.at(i);
        i_names.push_back(ss.str());
      }
    }

    // Adjoint seeds
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<nOut(); ++i) {
        ss.str(string());
        ss << "adj" << d << "_" << oscheme.at(i);
        i_names.push_back(ss.str());
      }
    }

    // Names of outputs
    std::vector<std::string> o_names;
    o_names.reserve(nOut()*(1+nfwd)+nIn()*nadj);

    // Nondifferentiated inputs
    for (int i=0; i<nOut(); ++i) {
      o_names.push_back("der_" + oscheme.at(i));
    }

    // Forward sensitivities
    for (int d=0; d<nfwd; ++d) {
      for (int i=0; i<nOut(); ++i) {
        ss.str(string());
        ss << "fwd" << d << "_" << oscheme.at(i);
        o_names.push_back(ss.str());
      }
    }

    // Adjoint sensitivities
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<nIn(); ++i) {
        ss.str(string());
        ss << "adj" << d << "_" << ischeme.at(i);
        o_names.push_back(ss.str());
      }
    }

    // Construct return function
    MXFunction ret(ss.str(), ret_in, ret_out,
                   make_dict("input_scheme", i_names, "output_scheme", o_names));

    // Consistency check for inputs
    int ind=0;
    for (int d=-1; d<nfwd; ++d) {
      for (int i=0; i<nIn(); ++i, ++ind) {
        if (ret.input(ind).nnz()!=0 && ret.input(ind).sparsity()!=input(i).sparsity()) {
          casadi_error("Incorrect sparsity for " << ret << " input " << ind << " \""
                       << i_names.at(ind) << "\". Expected " << input(i).dimString()
                       << " but got " << ret.input(ind).dimString());
        }
      }
    }
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<nOut(); ++i, ++ind) {
        if (ret.input(ind).nnz()!=0 && ret.input(ind).sparsity()!=output(i).sparsity()) {
          casadi_error("Incorrect sparsity for " << ret << " input " << ind <<
                       " \"" << i_names.at(ind) << "\". Expected " << output(i).dimString()
                       << " but got " << ret.input(ind).dimString());
        }
      }
    }

    // Consistency check for outputs
    ind=0;
    for (int d=-1; d<nfwd; ++d) {
      for (int i=0; i<nOut(); ++i, ++ind) {
        if (ret.output(ind).nnz()!=0 && ret.output(ind).sparsity()!=output(i).sparsity()) {
          casadi_error("Incorrect sparsity for " << ret << " output " << ind <<
                       " \"" <<  o_names.at(ind) << "\". Expected " << output(i).dimString()
                       << " but got " << ret.output(ind).dimString());
        }
      }
    }
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<nIn(); ++i, ++ind) {
        if (ret.output(ind).nnz()!=0 && ret.output(ind).sparsity()!=input(i).sparsity()) {
          casadi_error("Incorrect sparsity for " << ret << " output " << ind << " \""
                       << o_names.at(ind) << "\". Expected " << input(i).dimString()
                       << " but got " << ret.output(ind).dimString());
        }
      }
    }
    return ret;
  }

  Function Function::derForward(int nfwd) {
    return (*this)->derForward(nfwd);
  }

  Function Function::derReverse(int nadj) {
    return (*this)->derReverse(nadj);
  }

  void Function::setDerForward(const Function& fcn, int nfwd) {
    (*this)->setDerForward(fcn, nfwd);
  }

  void Function::setDerReverse(const Function& fcn, int nadj) {
    (*this)->setDerReverse(fcn, nadj);
  }

  void Function::printDimensions(std::ostream &stream) const {
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

  void Function::callDerivative(const DMatrixVector& arg, DMatrixVector& res,
                          const DMatrixVectorVector& fseed, DMatrixVectorVector& fsens,
                          const DMatrixVectorVector& aseed, DMatrixVectorVector& asens,
                          bool always_inline, bool never_inline) {
    call(arg, res, always_inline, never_inline);
    callForward(arg, res, fseed, fsens, always_inline, never_inline);
    callReverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  void Function::callDerivative(const SXVector& arg, SXVector& res,
                          const SXVectorVector& fseed, SXVectorVector& fsens,
                          const SXVectorVector& aseed, SXVectorVector& asens,
                          bool always_inline, bool never_inline) {
    call(arg, res, always_inline, never_inline);
    callForward(arg, res, fseed, fsens, always_inline, never_inline);
    callReverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  void Function::callDerivative(const MXVector& arg, MXVector& res,
                          const MXVectorVector& fseed, MXVectorVector& fsens,
                          const MXVectorVector& aseed, MXVectorVector& asens,
                          bool always_inline, bool never_inline) {
    call(arg, res, always_inline, never_inline);
    callForward(arg, res, fseed, fsens, always_inline, never_inline);
    callReverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  std::string Function::getSanitizedName() const {
    return (*this)->getSanitizedName();
  }

  std::string Function::sanitizeName(const std::string& name) {
    return FunctionInternal::sanitizeName(name);
  }

  void Function::callForward(const std::vector<MX>& arg, const std::vector<MX>& res,
                         const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens,
                         bool always_inline, bool never_inline) {
    checkArg(arg);
    checkRes(res);
    if (!matchingFwdSeed(fseed)) {
      return callForward(arg, res, replaceFwdSeed(fseed), fsens,
                     always_inline, never_inline);
    }
    (*this)->callForward(arg, res, fseed, fsens, always_inline, never_inline);
  }

  void Function::callReverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                         const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens,
                         bool always_inline, bool never_inline) {
    checkArg(arg);
    checkRes(res);
    if (!matchingAdjSeed(aseed)) {
      return callReverse(arg, res, replaceAdjSeed(aseed), asens,
                     always_inline, never_inline);
    }
    (*this)->callReverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  void Function::callForward(const std::vector<SX>& arg, const std::vector<SX>& res,
                         const std::vector<std::vector<SX> >& fseed,
                         std::vector<std::vector<SX> >& fsens,
                         bool always_inline, bool never_inline) {
    checkArg(arg);
    checkRes(res);
    if (!matchingFwdSeed(fseed)) {
      return callForward(arg, res, replaceFwdSeed(fseed), fsens,
                     always_inline, never_inline);
    }
    (*this)->callForward(arg, res, fseed, fsens, always_inline, never_inline);
  }

  void Function::callReverse(const std::vector<SX>& arg, const std::vector<SX>& res,
                         const std::vector<std::vector<SX> >& aseed,
                         std::vector<std::vector<SX> >& asens,
                         bool always_inline, bool never_inline) {
    checkArg(arg);
    checkRes(res);
    if (!matchingAdjSeed(aseed)) {
      return callReverse(arg, res, replaceAdjSeed(aseed), asens,
                     always_inline, never_inline);
    }
    (*this)->callReverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  void Function::callForward(const std::vector<DMatrix>& arg, const std::vector<DMatrix>& res,
                         const std::vector<std::vector<DMatrix> >& fseed,
                         std::vector<std::vector<DMatrix> >& fsens,
                         bool always_inline, bool never_inline) {
    checkArg(arg);
    checkRes(res);
    if (!matchingFwdSeed(fseed)) {
      return callForward(arg, res, replaceFwdSeed(fseed), fsens,
                     always_inline, never_inline);
    }
    (*this)->callForward(arg, res, fseed, fsens, always_inline, never_inline);
  }

  void Function::callReverse(const std::vector<DMatrix>& arg, const std::vector<DMatrix>& res,
                         const std::vector<std::vector<DMatrix> >& aseed,
                         std::vector<std::vector<DMatrix> >& asens,
                         bool always_inline, bool never_inline) {
    checkArg(arg);
    checkRes(res);
    if (!matchingAdjSeed(aseed)) {
      return callReverse(arg, res, replaceAdjSeed(aseed), asens,
                     always_inline, never_inline);
    }
    (*this)->callReverse(arg, res, aseed, asens, always_inline, never_inline);
  }

  std::vector<DMatrix> Function::operator()(const std::vector<DMatrix>& arg,
                                            bool always_inline, bool never_inline) {
    std::vector<DMatrix> res;
    call(arg, res, always_inline, never_inline);
    return res;
  }

  std::vector<SX> Function::operator()(const std::vector<SX>& arg,
                                       bool always_inline, bool never_inline) {
    std::vector<SX> res;
    call(arg, res, always_inline, never_inline);
    return res;
  }

  std::vector<MX> Function::operator()(const std::vector<MX>& arg,
                                       bool always_inline, bool never_inline) {
    std::vector<MX> res;
    call(arg, res, always_inline, never_inline);
    return res;
  }

  template<typename M>
  const std::map<std::string, M>
  Function::callMap(const std::map<std::string, M>& arg, bool always_inline, bool never_inline) {
    casadi_assert(isInit());

    // Get default inputs
    std::vector<M> v(nIn());
    for (int i=0; i<v.size(); ++i) {
      v[i] = defaultInput(i);
    }

    // Assign provided inputs
    for (typename std::map<std::string, M>::const_iterator i=arg.begin(); i!=arg.end(); ++i) {
      v.at(inputIndex(i->first)) = i->second;
    }

    // Make call
    v = (*this)(v, always_inline, never_inline);

    // Save to map
    std::map<std::string, M> ret;
    for (int i=0; i<v.size(); ++i) {
      ret[outputName(i)] = v[i];
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
    return arg.shape()==inp.shape() || arg.isempty() || arg.isscalar() ||
      (inp.size2()==arg.size1() && inp.size1()==arg.size2()
       && (arg.iscolumn() || inp.iscolumn())) ||
      (hcat && arg.size1()==inp.size1() && arg.size2() % inp.size2()==0);
  }

  template<typename M>
  void Function::checkArg(const std::vector<M>& arg, bool hcat) const {
    int n_in = nIn();
    casadi_assert_message(arg.size()==n_in, "Incorrect number of inputs: Expected "
                          << n_in << ", got " << arg.size());
    for (int i=0; i<n_in; ++i) {
      casadi_assert_message(checkMat(arg[i].sparsity(), input(i).sparsity(), hcat),
                            "Input " << i << " has mismatching shape. Expected "
                            << input(i).shape() << ", got " << arg[i].shape());
    }
  }

  template<typename M>
  void Function::checkRes(const std::vector<M>& res) const {
    int n_out = nOut();
    casadi_assert_message(res.size()==n_out, "Incorrect number of outputs: Expected "
                          << n_out << ", got " << res.size());
    for (int i=0; i<n_out; ++i) {
      casadi_assert_message(checkMat(res[i].sparsity(), output(i).sparsity()),
                            "Output " << i << " has mismatching shape. Expected "
                            << output(i).shape() << ", got " << res[i].shape());
    }
  }

  template<typename M>
  void Function::checkFwdSeed(const std::vector<std::vector<M> >& fseed) const {
    int n_in = nIn();
    for (int d=0; d<fseed.size(); ++d) {
      casadi_assert_message(fseed[d].size()==n_in,
                            "Incorrect number of forward seeds for direction " << d
                            << ": Expected " << n_in << ", got " << fseed[d].size());
      for (int i=0; i<n_in; ++i) {
        casadi_assert_message(checkMat(fseed[d][i].sparsity(), input(i).sparsity()),
                              "Forward seed " << i << " for direction " << d
                              << " has mismatching shape. Expected " << input(i).shape()
                              << ", got " << fseed[d][i].shape());
      }
    }
  }

  template<typename M>
  void Function::checkAdjSeed(const std::vector<std::vector<M> >& aseed) const {
    int n_out = nOut();
    for (int d=0; d<aseed.size(); ++d) {
      casadi_assert_message(aseed[d].size()==n_out,
                            "Incorrect number of adjoint seeds for direction " << d
                            << ": Expected " << n_out << ", got " << aseed[d].size());
      for (int i=0; i<n_out; ++i) {
        casadi_assert_message(checkMat(aseed[d][i].sparsity(), output(i).sparsity()),
                              "Adjoint seed " << i << " for direction " << d
                              << " has mismatching shape. Expected " << output(i).shape()
                              << ", got " << aseed[d][i].shape());
      }
    }
  }

  template<typename M>
  bool Function::matchingArg(const std::vector<M>& arg, bool hcat) const {
    checkArg(arg, hcat);
    int n_in = nIn();
    for (int i=0; i<n_in; ++i) {
      if (hcat) {
        if (arg.at(i).size1()!=input(i).size1()) return false;
        if (arg.at(i).size2() % input(i).size2()!=0 || arg.at(i).size2()==0) return false;
      } else {
        if (arg.at(i).shape()!=input(i).shape()) return false;
      }
    }
    return true;
  }

  template<typename M>
  bool Function::matchingRes(const std::vector<M>& res) const {
    checkRes(res);
    int n_out = nOut();
    for (int i=0; i<n_out; ++i) {
      if (res.at(i).shape()!=output(i).shape()) return false;
    }
    return true;
  }

  template<typename M>
  bool Function::matchingFwdSeed(const std::vector<std::vector<M> >& fseed) const {
    checkFwdSeed(fseed);
    for (int d=0; d<fseed.size(); ++d) {
      if (!matchingArg(fseed[d])) return false;
    }
    return true;
  }

  template<typename M>
  bool Function::matchingAdjSeed(const std::vector<std::vector<M> >& aseed) const {
    checkAdjSeed(aseed);
    for (int d=0; d<aseed.size(); ++d) {
      if (!matchingRes(aseed[d])) return false;
    }
    return true;
  }

  template<typename M>
  M replaceMat(const M& arg, const Sparsity& inp, bool hcat=false) {
    if (arg.shape()==inp.shape()) {
      // Matching dimensions already
      return arg;
    } else if (hcat && arg.size1()==inp.size1() && arg.size2() % inp.size2()==0
                    && arg.size2() >=0) {
      // Matching horzcat dimensions
      return arg;
    } else if (arg.isempty()) {
      // Empty matrix means set zero
      return M(inp.shape());
    } else if (arg.isscalar()) {
      // Scalar assign means set all
      return M(inp, arg);
    } else {
      // Assign vector with transposing
      casadi_assert(arg.size1()==inp.size2() && arg.size2()==inp.size1()
                    && (arg.iscolumn() || inp.iscolumn()));
      return arg.T();
    }
  }

  template<typename M>
  std::vector<M> Function::replaceArg(const std::vector<M>& arg, bool hcat) const {
    std::vector<M> r(arg.size());
    for (int i=0; i<r.size(); ++i) r[i] = replaceMat(arg[i], input(i).sparsity(), hcat);
    return r;
  }

  template<typename M>
  std::vector<M> Function::replaceRes(const std::vector<M>& res) const {
    std::vector<M> r(res.size());
    for (int i=0; i<r.size(); ++i) r[i] = replaceMat(res[i], output(i).sparsity());
    return r;
  }

  template<typename M>
  std::vector<std::vector<M> >
  Function::replaceFwdSeed(const std::vector<std::vector<M> >& fseed) const {
    std::vector<std::vector<M> > r(fseed.size());
    for (int d=0; d<r.size(); ++d) r[d] = replaceArg(fseed[d]);
    return r;
  }

  template<typename M>
  std::vector<std::vector<M> >
  Function::replaceAdjSeed(const std::vector<std::vector<M> >& aseed) const {
    std::vector<std::vector<M> > r(aseed.size());
    for (int d=0; d<r.size(); ++d) r[d] = replaceRes(aseed[d]);
    return r;
  }

  double Function::defaultInput(int ind) const {
    return (*this)->defaultInput(ind);
  }

  void Function::operator()(const double** arg, double** res, int* iw, double* w) {
    (*this)->evalD(arg, res, iw, w);
  }

} // namespace casadi

