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
#include "../mx/call_function.hpp"
#include "../function/mx_function.hpp"
#include <typeinfo>
#include "../std_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "parallelizer.hpp"

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

  vector<DMatrix> Function::call(const vector<DMatrix> &arg,
                                 bool always_inline, bool never_inline) {
    DMatrixVectorVector dummy;
    DMatrixVector res;
    callDerivative(arg, res, dummy, dummy, dummy, dummy, always_inline, never_inline);
    return res;
  }

  vector<SX> Function::call(const vector<SX> &arg, bool always_inline, bool never_inline) {
    SXVectorVector dummy;
    SXVector res;
    callDerivative(arg, res, dummy, dummy, dummy, dummy, always_inline, never_inline);
    return res;
  }

  vector<MX> Function::call(const vector<MX> &arg, bool always_inline, bool never_inline) {
    MXVectorVector dummy;
    MXVector res;
    callDerivative(arg, res, dummy, dummy, dummy, dummy, always_inline, never_inline);
    return res;
  }

  vector<vector<MX> > Function::callParallel(const vector<vector<MX> > &x,
                                             const Dictionary& paropt) {
    assertInit();

    // Make sure not empty
    casadi_assert_message(x.size()>1, "Function: callParallel(vector<vector<MX> >): "
                          "argument must be of length > 1. You supplied length "
                          << x.size() << ".");

    // Return object
    vector<vector<MX> > ret(x.size());

    // Check if we are bypassing the parallelizer
    Dictionary::const_iterator ii=paropt.find("parallelization");
    if (ii!=paropt.end() && ii->second=="expand") {
      for (int i=0; i<x.size(); ++i) {
        ret[i] = call(x[i]);
      }
      return ret;
    }

    // Create parallelizer object and initialize it
    Parallelizer p(vector<Function>(x.size(), *this));
    p.setOption(paropt);
    p.init();

    // Concatenate the arguments
    vector<MX> p_in;
    p_in.reserve(x.size() * getNumInputs());
    for (int i=0; i<x.size(); ++i) {
      p_in.insert(p_in.end(), x[i].begin(), x[i].end());
      p_in.resize(p_in.size()+getNumInputs()-x[i].size());
    }

    // Call the parallelizer
    vector<MX> p_out = p.call(p_in);
    casadi_assert(p_out.size() == x.size() * getNumOutputs());

    // Collect the outputs
    vector<MX>::const_iterator it=p_out.begin();
    for (int i=0; i<x.size(); ++i) {
      ret[i].insert(ret[i].end(), it, it+getNumOutputs());
      it += getNumOutputs();
    }
    return ret;
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

  bool Function::spCanEvaluate(bool fwd) {
    return (*this)->spCanEvaluate(fwd);
  }

  void Function::spInit(bool fwd) {
    (*this)->spInit(fwd);
  }

  Function Function::derivative(int nfwd, int nadj) {
    return (*this)->derivative(nfwd, nadj);
  }

  void Function::setDerivative(const Function& fcn, int nfwd, int nadj) {
    (*this)->setDerivative(fcn, nfwd, nadj);
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

  std::string Function::generateCode() {
    std::ostringstream cfile;
    generateCode(cfile);
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
    (*this)->call(arg, res, fseed, fsens, aseed, asens, always_inline, never_inline);
  }

  void Function::callDerivative(const SXVector& arg, SXVector& res,
                          const SXVectorVector& fseed, SXVectorVector& fsens,
                          const SXVectorVector& aseed, SXVectorVector& asens,
                          bool always_inline, bool never_inline) {
    casadi_assert_message(arg.size()==getNumInputs(), "Function::callDerivative: dimension "
                          "mismatch. You supplied " << arg.size()
                          << " arguments instead of expected " << getNumInputs() << ".");
    (*this)->call(arg, res, fseed, fsens, aseed, asens, always_inline, never_inline);
  }

  void Function::callDerivative(const MXVector& arg, MXVector& res,
                          const MXVectorVector& fseed, MXVectorVector& fsens,
                          const MXVectorVector& aseed, MXVectorVector& asens,
                          bool always_inline, bool never_inline) {
    casadi_assert_message(arg.size()==getNumInputs(), "Function::callDerivative: "
                          "dimension mismatch. You supplied "
                          << arg.size() << " arguments instead of expected "
                          << getNumInputs() << ".");
    (*this)->call(arg, res, fseed, fsens, aseed, asens, always_inline, never_inline);
  }

  std::string Function::getSanitizedName() const {
    return (*this)->getSanitizedName();
  }
} // namespace casadi

