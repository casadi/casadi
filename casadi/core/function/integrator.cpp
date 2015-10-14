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


#include "integrator.hpp"
#include "integrator_internal.hpp"
#include "sx_function.hpp"
#include "mx_function.hpp"

using namespace std;
namespace casadi {

  Integrator::Integrator() {
  }

  Integrator::Integrator(const std::string& name, const std::string& solver, const Function& f,
                         const Dict& opts) {
    // Backwards DAE
    Function g;
    Dict opts2 = opts;
    Dict::const_iterator it=opts2.find("rdae");
    if (it!=opts2.end()) {
      g = it->second;
      opts2.erase(it);
    }

    // Create an initialize
    assignNode(IntegratorInternal::getPlugin(solver).creator(name, f, g));
    setOption(opts2);
    init();
  }

  Integrator::Integrator(const std::string& name, const std::string& solver,
                         const std::pair<Function, Function>& fg, const Dict& opts) {
    assignNode(IntegratorInternal::getPlugin(solver).creator(name, fg.first, fg.second));
    setOption(opts);
    init();
  }

  Integrator::Integrator(const std::string& name, const std::string& solver,
                         const SXDict& dae, const Dict& opts) {
    SX x, z, p, t, ode, alg, quad;
    for (SXDict::const_iterator i=dae.begin(); i!=dae.end(); ++i) {
      if (i->first=="x") {
        x = i->second;
      } else if (i->first=="z") {
        z = i->second;
      } else if (i->first=="p") {
        p = i->second;
      } else if (i->first=="t") {
        t = i->second;
      } else if (i->first=="ode") {
        ode = i->second;
      } else if (i->first=="alg") {
        alg = i->second;
      } else if (i->first=="quad") {
        quad = i->second;
      } else {
        casadi_error("No such field: \"" + i->first + "\"");
      }
    }
    Function f=SX::fun("dae", daeIn("x", x, "z", z, "p", p),
                    daeOut("ode", ode, "alg", alg, "quad", quad));
    Function g; // dummy

    // Load the plugin
    assignNode(IntegratorInternal::getPlugin(solver).creator(name, f, g));
    setOption(opts);
    init();
  }

  Integrator::Integrator(const std::string& name, const std::string& solver,
                         const MXDict& dae, const Dict& opts) {
    MX x, z, p, t, ode, alg, quad;
    for (MXDict::const_iterator i=dae.begin(); i!=dae.end(); ++i) {
      if (i->first=="x") {
        x = i->second;
      } else if (i->first=="z") {
        z = i->second;
      } else if (i->first=="p") {
        p = i->second;
      } else if (i->first=="t") {
        t = i->second;
      } else if (i->first=="ode") {
        ode = i->second;
      } else if (i->first=="alg") {
        alg = i->second;
      } else if (i->first=="quad") {
        quad = i->second;
      } else {
        casadi_error("No such field: \"" + i->first + "\"");
      }
    }
    Function f=MX::fun("dae", daeIn("x", x, "z", z, "p", p),
                    daeOut("ode", ode, "alg", alg, "quad", quad));
    Function g; // dummy

    // Load the plugin
    assignNode(IntegratorInternal::getPlugin(solver).creator(name, f, g));
    setOption(opts);
    init();
  }

  void Integrator::printStats(ostream &stream) const {
    (*this)->printStats(stream);
  }

  IntegratorInternal* Integrator::operator->() {
    return static_cast<IntegratorInternal*>(Function::operator->());
  }

  const IntegratorInternal* Integrator::operator->() const {
    return static_cast<const IntegratorInternal*>(Function::operator->());
  }

  void Integrator::reset() {
    (*this)->reset();
  }

  void Integrator::integrate(double t_out) {
    (*this)->integrate(t_out);
  }

  bool Integrator::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const IntegratorInternal*>(ptr)!=0;
  }

  void Integrator::resetB() {
    (*this)->resetB();
  }

  void Integrator::integrateB(double t_out) {
    (*this)->integrateB(t_out);
  }

  Function Integrator::getDAE() {
    return (*this)->f_;
  }

  std::pair<Function, Function> Integrator::getAugmented(int nfwd, int nadj) {
    IntegratorInternal::AugOffset offset;
    return (*this)->getAugmented(nfwd, nadj, offset);
  }

  bool Integrator::hasPlugin(const std::string& name) {
    return IntegratorInternal::hasPlugin(name);
  }

  void Integrator::loadPlugin(const std::string& name) {
    IntegratorInternal::loadPlugin(name);
  }

  std::string Integrator::doc(const std::string& name) {
    return IntegratorInternal::getPlugin(name).doc;
  }

  void Integrator::setStopTime(double tf) {
    (*this)->setStopTime(tf);
  }

} // namespace casadi

