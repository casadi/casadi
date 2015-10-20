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


#include "integrator_internal.hpp"

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

    Integrator tmp(name, solver, make_pair(f, g), opts2);
    assignNode(tmp.get());
  }

  Integrator::Integrator(const std::string& name, const std::string& solver,
                         const std::pair<Function, Function>& fg, const Dict& opts) {
    Function f=fg.first, g=fg.second;
    if (f.is_a("sxfunction")) {
      SXDict dae;
      vector<SX> v = SX::get_input(f), vf=v, vg=v;
      for (auto&& i : v) casadi_assert(i.isSymbolic());
      dae["t"] = v[DAE_T];
      dae["x"] = v[DAE_X];
      dae["z"] = v[DAE_Z];
      dae["p"] = v[DAE_P];
      v = f(v);
      dae["ode"] = v[DAE_ODE];
      dae["alg"] = v[DAE_ALG];
      dae["quad"] = v[DAE_QUAD];
      if (!g.isNull()) {
        v = SX::get_input(g);
        for (auto&& i : v) casadi_assert(i.isSymbolic());
        dae["rx"] = v[RDAE_RX];
        dae["rz"] = v[RDAE_RZ];
        dae["rp"] = v[RDAE_RP];
        vg[DAE_T] = v[RDAE_T];
        vg[DAE_X] = v[RDAE_X];
        vg[DAE_Z] = v[RDAE_Z];
        vg[DAE_P] = v[RDAE_P];
        v = substitute(g(v), vg, vf);
        dae["rode"] = v[RDAE_ODE];
        dae["ralg"] = v[RDAE_ALG];
        dae["rquad"] = v[RDAE_QUAD];
      }
      Function tmp = Function::integrator(name, solver, dae, opts);
      assignNode(tmp.get());
    } else {
      casadi_assert(f.is_a("mxfunction"));
      MXDict dae;
      vector<MX> v = MX::get_input(f), vf=v, vg=v;
      for (auto&& i : v) casadi_assert(i.isSymbolic());
      dae["t"] = v[DAE_T];
      dae["x"] = v[DAE_X];
      dae["z"] = v[DAE_Z];
      dae["p"] = v[DAE_P];
      v = f(v);
      dae["ode"] = v[DAE_ODE];
      dae["alg"] = v[DAE_ALG];
      dae["quad"] = v[DAE_QUAD];
      if (!g.isNull()) {
        v = MX::get_input(g);
        for (auto&& i : v) casadi_assert(i.isSymbolic());
        dae["rx"] = v[RDAE_RX];
        dae["rz"] = v[RDAE_RZ];
        dae["rp"] = v[RDAE_RP];
        vg[DAE_T] = v[RDAE_T];
        vg[DAE_X] = v[RDAE_X];
        vg[DAE_Z] = v[RDAE_Z];
        vg[DAE_P] = v[RDAE_P];
        v = substitute(g(v), vg, vf);
        dae["rode"] = v[RDAE_ODE];
        dae["ralg"] = v[RDAE_ALG];
        dae["rquad"] = v[RDAE_QUAD];
      }
      Function tmp = Function::integrator(name, solver, dae, opts);
      assignNode(tmp.get());
    }
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

  bool Integrator::test_cast(const SharedObjectNode* ptr) {
    return dynamic_cast<const IntegratorInternal*>(ptr)!=0;
  }

  void Integrator::resetB() {
    (*this)->resetB();
  }

  void Integrator::integrateB(double t_out) {
    (*this)->integrateB(t_out);
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

