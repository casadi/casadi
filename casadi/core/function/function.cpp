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
#include "../matrix/matrix_tools.hpp"

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

    // Call the internal function
    return (*this)->createMap(x, parallelization);
  }

  void Function::evaluate() {
    assertInit();
    (*this)->evaluate();
  }

  int Function::getNumInputNonzeros() const {
    return (*this)->getNumInputNonzeros();
  }

  int Function::getNumOutputNonzeros() const {
    return (*this)->getNumOutputNonzeros();
  }

  int Function::getNumInputElements() const {
    return (*this)->getNumInputElements();
  }

  int Function::getNumOutputElements() const {
    return (*this)->getNumOutputElements();
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

  const Dictionary & Function::getStats() const {
    return (*this)->getStats();
  }

  GenericType Function::getStat(const string& name) const {
    return (*this)->getStat(name);
  }

  Sparsity& Function::jacSparsity(int iind, int oind, bool compact, bool symmetric) {
    return (*this)->jacSparsity(iind, oind, compact, symmetric);
  }

  void Function::setJacSparsity(const Sparsity& sp, int iind, int oind, bool compact) {
    (*this)->setJacSparsity(sp, iind, oind, compact);
  }

  std::vector<MX> Function::symbolicInput() const {
    return (*this)->symbolicInput();
  }

  std::vector<SX> Function::symbolicInputSX() const {
    return (*this)->symbolicInputSX();
  }

  void Function::setInputScheme(const IOScheme &scheme) {
    return (*this)->setInputScheme(scheme);
  }

  void Function::setOutputScheme(const IOScheme &scheme) {
    return (*this)->setOutputScheme(scheme);
  }

  IOScheme Function::getInputScheme() const {
    return (*this)->getInputScheme();
  }

  IOScheme Function::getOutputScheme() const {
    return (*this)->getOutputScheme();
  }

  void Function::spEvaluate(bool fwd) {
    (*this)->spEvaluate(fwd);
  }

  void Function::nTmp(size_t& ni, size_t& nr) {
    (*this)->nTmp(ni, nr);
  }

  void Function::spFwd(cp_bvec_t* arg, p_bvec_t* res, int* itmp, bvec_t* rtmp) {
    (*this)->spFwdSwitch(arg, res, itmp, rtmp);
  }

  void Function::spAdj(p_bvec_t* arg, p_bvec_t* res, int* itmp, bvec_t* rtmp) {
    (*this)->spAdjSwitch(arg, res, itmp, rtmp);
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
    int num_in = getNumInputs();
    int num_out = getNumOutputs();

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
          *it = it->setSparse(output(i).sparsity());
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
          *it = it->setSparse(input(i).sparsity());
      ret_out.insert(ret_out.end(), res.begin(), res.end());
    }

    // Construct return function
    MXFunction ret(ret_in, ret_out);

    // Give it a suitable name
    stringstream ss;
    ss << "derivative_" << getOption("name") << "_" << nfwd << "_" << nadj;
    ret.setOption("name", ss.str());

    // Names of inputs
    std::vector<std::string> i_names;
    i_names.reserve(getNumInputs()*(1+nfwd)+getNumOutputs()*nadj);

    // Nondifferentiated inputs
    for (int i=0; i<getNumInputs(); ++i) {
      i_names.push_back("der_" + inputScheme().entryLabel(i));
    }

    // Forward seeds
    for (int d=0; d<nfwd; ++d) {
      for (int i=0; i<getNumInputs(); ++i) {
        ss.str(string());
        ss << "fwd" << d << "_" << inputScheme().entryLabel(i);
        i_names.push_back(ss.str());
      }
    }

    // Adjoint seeds
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<getNumOutputs(); ++i) {
        ss.str(string());
        ss << "adj" << d << "_" << outputScheme().entryLabel(i);
        i_names.push_back(ss.str());
      }
    }

    // Pass to return object
    ret.setInputScheme(i_names);

    // Names of outputs
    std::vector<std::string> o_names;
    o_names.reserve(getNumOutputs()*(1+nfwd)+getNumInputs()*nadj);

    // Nondifferentiated inputs
    for (int i=0; i<getNumOutputs(); ++i) {
      o_names.push_back("der_" + outputScheme().entryLabel(i));
    }

    // Forward sensitivities
    for (int d=0; d<nfwd; ++d) {
      for (int i=0; i<getNumOutputs(); ++i) {
        ss.str(string());
        ss << "fwd" << d << "_" << outputScheme().entryLabel(i);
        o_names.push_back(ss.str());
      }
    }

    // Adjoint sensitivities
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<getNumInputs(); ++i) {
        ss.str(string());
        ss << "adj" << d << "_" << inputScheme().entryLabel(i);
        o_names.push_back(ss.str());
      }
    }

    // Pass to return object
    ret.setOutputScheme(o_names);

    // Initialize it
    ret.init();

    // Consistency check for inputs
    int ind=0;
    for (int d=-1; d<nfwd; ++d) {
      for (int i=0; i<getNumInputs(); ++i, ++ind) {
        if (ret.input(ind).nnz()!=0 && ret.input(ind).sparsity()!=input(i).sparsity()) {
          casadi_error("Incorrect sparsity for " << ret << " input " << ind << " \""
                       << i_names.at(ind) << "\". Expected " << input(i).dimString()
                       << " but got " << ret.input(ind).dimString());
        }
      }
    }
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<getNumOutputs(); ++i, ++ind) {
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
      for (int i=0; i<getNumOutputs(); ++i, ++ind) {
        if (ret.output(ind).nnz()!=0 && ret.output(ind).sparsity()!=output(i).sparsity()) {
          casadi_error("Incorrect sparsity for " << ret << " output " << ind <<
                       " \"" <<  o_names.at(ind) << "\". Expected " << output(i).dimString()
                       << " but got " << ret.output(ind).dimString());
        }
      }
    }
    for (int d=0; d<nadj; ++d) {
      for (int i=0; i<getNumInputs(); ++i, ++ind) {
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

  void Function::generateCode(const string& filename, bool generate_main) {
    // Detect a C++ ending
    vector<string> cpp_endings;
    cpp_endings.push_back(".cpp");
    cpp_endings.push_back(".cxx");
    cpp_endings.push_back(".cc");
    cpp_endings.push_back(".cp");
    cpp_endings.push_back(".c++");
    for (vector<string>::const_iterator it=cpp_endings.begin(); it!=cpp_endings.end(); ++it) {
      if (filename.size()>it->size() &&
         filename.compare(filename.size()-it->size(), it->size(), *it)==0) {
        casadi_warning("Function::generateCode: Detected C++ file ending "
                       "(generated code is C, not C++)");
      }
    }

    // Create a file
    std::ofstream cfile;
    cfile.open(filename.c_str());
    generateCode(cfile, generate_main);
    cfile.close();
  }

  void Function::generateFunction(std::ostream &stream, const std::string& fname,
                                  const std::string& type, CodeGenerator& gen) const {
    (*this)->generateFunction(stream, fname, type, gen);
  }

  std::string Function::generateCodeStr(bool generate_main) {
    std::ostringstream cfile;
    generateCode(cfile, generate_main);
    return cfile.str();
  }

  void Function::generateCode(std::ostream &stream, bool generate_main) {
    (*this)->generateCode(stream, generate_main);
  }

  const IOScheme& Function::inputScheme() const {
    return (*this)->inputScheme();
  }

  const IOScheme& Function::outputScheme() const {
    return (*this)->outputScheme();
  }

  IOScheme& Function::inputScheme() {
    return (*this)->inputScheme();
  }

  IOScheme& Function::outputScheme() {
    return (*this)->outputScheme();
  }

  const IOSchemeVector<DMatrix>& Function::input_struct() const {
    return (*this)->input_struct();
  }

  const IOSchemeVector<DMatrix>& Function::output_struct() const {
    return (*this)->output_struct();
  }

  IOSchemeVector<DMatrix>& Function::input_struct() {
    return (*this)->input_struct();
  }

  IOSchemeVector<DMatrix>& Function::output_struct() {
    return (*this)->output_struct();
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

  IOSchemeVector<DMatrix> Function::
  operator()(const IOSchemeVector<DMatrix>& arg, bool always_inline, bool never_inline) {
    return outputScheme().fromVector(operator()(arg.data, always_inline, never_inline));
  }

  IOSchemeVector<SX> Function::
  operator()(const IOSchemeVector<SX>& arg, bool always_inline, bool never_inline) {
    return outputScheme().fromVector(operator()(arg.data, always_inline, never_inline));
  }

  IOSchemeVector<MX> Function::
  operator()(const IOSchemeVector<MX>& arg, bool always_inline, bool never_inline) {
    return outputScheme().fromVector(operator()(arg.data, always_inline, never_inline));
  }

  inline bool checkMat(const Sparsity& arg, const Sparsity& inp) {
    return arg.shape()==inp.shape() || arg.isEmpty() || arg.isScalar() ||
      (inp.size2()==arg.size1() && inp.size1()==arg.size2()
       && (arg.isVector() || inp.isVector()));
  }

  template<typename M>
  void Function::checkArg(const std::vector<M>& arg) const {
    int n_in = getNumInputs();
    casadi_assert_message(arg.size()==n_in, "Incorrect number of inputs: Expected "
                          << n_in << ", got " << arg.size());
    for (int i=0; i<n_in; ++i) {
      casadi_assert_message(checkMat(arg[i].sparsity(), input(i).sparsity()),
                            "Input " << i << " has mismatching shape. Expected "
                            << input(i).shape() << ", got " << arg[i].shape());
    }
  }

  template<typename M>
  void Function::checkRes(const std::vector<M>& res) const {
    int n_out = getNumOutputs();
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
    int n_in = getNumInputs();
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
    int n_out = getNumOutputs();
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
  bool Function::matchingArg(const std::vector<M>& arg) const {
    checkArg(arg);
    int n_in = getNumInputs();
    for (int i=0; i<n_in; ++i) {
      if (arg.at(i).shape()!=input(i).shape()) return false;
    }
    return true;
  }

  template<typename M>
  bool Function::matchingRes(const std::vector<M>& res) const {
    checkRes(res);
    int n_out = getNumOutputs();
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
  M replaceMat(const M& arg, const Sparsity& inp) {
    if (arg.shape()==inp.shape()) {
      // Matching dimensions already
      return arg;
    } else if (arg.isEmpty()) {
      // Empty matrix means set zero
      return M(inp.shape());
    } else if (arg.isScalar()) {
      // Scalar assign means set all
      return M(inp, arg);
    } else {
      // Assign vector with transposing
      casadi_assert(arg.size1()==inp.size2() && arg.size2()==inp.size1()
                    && (arg.isVector() || inp.isVector()));
      return arg.T();
    }
  }

  template<typename M>
  std::vector<M> Function::replaceArg(const std::vector<M>& arg) const {
    std::vector<M> r(arg.size());
    for (int i=0; i<r.size(); ++i) r[i] = replaceMat(arg[i], input(i).sparsity());
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

} // namespace casadi

