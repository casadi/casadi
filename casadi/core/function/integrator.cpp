/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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
#include <cassert>

// For dynamic loading
#ifdef WITH_DL
#ifdef _WIN32 // also for 64-bit
#include <windows.h>
#else // _WIN32
#include <dlfcn.h>
#endif // _WIN32
#endif // WITH_DL

using namespace std;
namespace casadi {

  Integrator::Integrator() {
  }

  Integrator::Integrator(const std::string& name, const Function& f, const Function& g) {
    // Check if the solver has been loaded
    std::map<std::string, IntegratorPlugin>::iterator it=solvers_.find(name);

    // Load the solver if needed
    if (it==solvers_.end()) {
      loadPlugin(name);
      it=solvers_.find(name);
    }
    casadi_assert(it!=solvers_.end());
    assignNode(it->second.creator(f, g));
  }

  std::map<std::string, Integrator::IntegratorPlugin> Integrator::solvers_;

  void Integrator::registerPlugin(RegFcn regfcn) {
    // Create a temporary struct
    IntegratorPlugin plugin;
   
    // Set the fields
    int flag = regfcn(&plugin);
    casadi_assert(flag==0);

    // Check if the solver name is in use
    std::map<std::string, IntegratorPlugin>::iterator it=solvers_.find(plugin.name);
    casadi_assert_message(it==solvers_.end(), "Solver " << plugin.name << " is already in use");

    // Add to list of solvers
    solvers_[plugin.name] = plugin;
  }

  void Integrator::loadPlugin(const std::string& name) {
#ifndef WITH_DL
    casadi_error("WITH_DL option needed for dynamic loading");
#else // WITH_DL
    std::string lib = "libcasadi_integrator_" + name + ".dylib";

    // Load the dll
    void* handle;
    RegFcn reg;
    std::string regName = "casadi_register_integrator_" + name;
#ifdef _WIN32
    handle = LoadLibrary(TEXT(lib.c_str()));
    casadi_assert_message(handle!=0, "Integrator: Cannot open function: "
                        << lib << ". error code (WIN32): "<< GetLastError());

    reg = (RegFcn)GetProcAddress(handle_, TEXT(regName.c_str()));
    if (reg==0) throw CasadiException("Integrator: no \"" + regName + "\" found");
#else // _WIN32
  handle  = dlopen(lib.c_str(), RTLD_LAZY);
  casadi_assert_message(handle!=0, "Integrator: Cannot open function: "
                        << lib << ". error code: "<< dlerror());

  // reset error
  dlerror();

  // Load creator
  reg = (RegFcn)dlsym(handle, regName.c_str());
  if (dlerror()) throw CasadiException("Integrator: no \"" + regName + "\" found");
#endif // _WIN32

  // Register the plugin
  registerPlugin(reg);

#endif // WITH_DL
  }

  Integrator  Integrator::clone() const {
    Integrator ret;
    if (!isNull()) ret.assignNode((*this)->clone());
    return ret;
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

  bool Integrator::checkNode() const {
    return dynamic_cast<const IntegratorInternal*>(get())!=0;
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

} // namespace casadi

