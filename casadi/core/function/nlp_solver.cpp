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


#include "nlp_solver.hpp"
#include "nlp_solver_internal.hpp"
#include "sx_function.hpp"
#include "mx_function.hpp"

using namespace std;
namespace casadi {

  NlpSolver::NlpSolver() {
  }

  NlpSolver::NlpSolver(const std::string& name, const std::string& solver,
                       const Function& nlp, const Dict& opts) {
    assignNode(NlpSolverInternal::instantiatePlugin(name, solver, nlp));
    setOption(opts);
    init(false);
  }

  NlpSolver::NlpSolver(const std::string& name,
                       const std::string& solver,
                       const SXDict& nlp,
                       const Dict& opts) {
    // Create an NLP function
    SX x, p, f, g;
    for (SXDict::const_iterator i=nlp.begin(); i!=nlp.end(); ++i) {
      if (i->first=="x") {
        x = i->second;
      } else if (i->first=="p") {
        p = i->second;
      } else if (i->first=="f") {
        f = i->second;
      } else if (i->first=="g") {
        g = i->second;
      } else {
        casadi_error("No such field: \"" + i->first + "\"");
      }
    }
    SXFunction nlpf("nlp", nlpIn("x", x, "p", p), nlpOut("f", f, "g", g));

    // Create the solver instance
    assignNode(NlpSolverInternal::instantiatePlugin(name, solver, nlpf));
    setOption(opts);
    init(false);
  }

  NlpSolver::NlpSolver(const std::string& name,
                       const std::string& solver,
                       const MXDict& nlp,
                       const Dict& opts) {
    // Create an NLP function
    MX x, p, f, g;
    for (MXDict::const_iterator i=nlp.begin(); i!=nlp.end(); ++i) {
      if (i->first=="x") {
        x = i->second;
      } else if (i->first=="p") {
        p = i->second;
      } else if (i->first=="f") {
        f = i->second;
      } else if (i->first=="g") {
        g = i->second;
      } else {
        casadi_error("No such field: \"" + i->first + "\"");
      }
    }
    MXFunction nlpf("nlp", nlpIn("x", x, "p", p), nlpOut("f", f, "g", g));

    // Create the solver instance
    assignNode(NlpSolverInternal::instantiatePlugin(name, solver, nlpf));
    setOption(opts);
    init(false);
  }

  NlpSolverInternal* NlpSolver::operator->() {
    return static_cast<NlpSolverInternal*>(Function::operator->());
  }

  const NlpSolverInternal* NlpSolver::operator->() const {
    return static_cast<const NlpSolverInternal*>(Function::operator->());
  }

  bool NlpSolver::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const NlpSolverInternal*>(ptr)!=0;
  }

  void NlpSolver::reportConstraints(std::ostream &stream) {
    (*this)->reportConstraints();
  }

  Function NlpSolver::nlp() {
    return (*this)->nlp_;
  }

  Function NlpSolver::gradF() {
    return (*this)->gradF();
  }

  Function NlpSolver::jacG() {
    return (*this)->jacG();
  }

  Function NlpSolver::hessLag() {
    return (*this)->hessLag();
  }

  bool NlpSolver::hasPlugin(const std::string& name) {
    return NlpSolverInternal::hasPlugin(name);
  }

  void NlpSolver::loadPlugin(const std::string& name) {
    NlpSolverInternal::loadPlugin(name);
  }

  std::string NlpSolver::doc(const std::string& name) {
    return NlpSolverInternal::getPlugin(name).doc;
  }

  DMatrix NlpSolver::getReducedHessian() {
    return (*this)->getReducedHessian();
  }

  void NlpSolver::setOptionsFromFile(const std::string & file) {
    (*this)->setOptionsFromFile(file);
  }

} // namespace casadi
