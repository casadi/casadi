/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl,
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


#include "erk.hpp"

namespace casadi {

  extern "C"
  int CASADI_INTEGRATOR_ERK_EXPORT
      casadi_register_integrator_erk(Integrator::Plugin* plugin) {
    plugin->creator = Erk::creator;
    plugin->name = "erk";
    plugin->doc = Erk::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Erk::options_;
    plugin->deserialize = &Erk::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_INTEGRATOR_ERK_EXPORT casadi_load_integrator_erk() {
    Integrator::registerPlugin(casadi_register_integrator_erk);
  }

  Erk::Erk(const std::string& name, const Function& dae, double t0,
      const std::vector<double>& tout) : RuntimeIntegrator(name, dae, t0, tout) {
  }

  Erk::~Erk() {
    clear_mem();
  }

  const Options Erk::options_
  = {{&RuntimeIntegrator::options_},
     {{"method",
       {OT_STRING,
        "Runge-Kutta method: 'rk45' (pair of order 5(4), adaptive unless step_size is set) or "
        "'rk4' (classical, order 4, fixed steps only; default step_size (tf - t0)/20) [rk45]"}}
     }
  };

  void Erk::init(const Dict& opts) {
    // Call the base class init
    RuntimeIntegrator::init(opts);

    // Default options
    method_ = "rk45";

    // Read options
    for (auto&& op : opts) {
      if (op.first == "method") {
        method_ = op.second.to_string();
      }
    }
    casadi_assert(nz_ == 0, "'erk' does not support algebraic states, use 'rosenbrock'");

    // Reverse communication problem
    set_coefficients();
    if (m_err_.empty() && step_size_ == 0) step_size_ = (tout_.back() - t0_) / 20;
    set_ivp_prob();
    set_nominal_prob();
    set_erk_prob();

    // Work vectors of the kernel
    alloc_ivp_work();
  }

  void Erk::set_coefficients() {
    if (method_ == "rk4") {
      // Classical method of order 4
      a_ = {0, 0, 0, 0,
            0.5, 0, 0, 0,
            0, 0.5, 0, 0,
            0, 0, 1, 0};
      alpha_ = {0, 0.5, 0.5, 1};
      m_ = {1./6, 1./3, 1./3, 1./6};
      // Dense output of order 3: m_i(theta), coefficients of theta, theta^2, theta^3
      n_dense_ = 3;
      m_dense_ = {1, 0, 0, 0,
                  -1.5, 1, 1, -0.5,
                  2./3, -2./3, -2./3, 2./3};
      // No error estimate
      m_err_.clear();
      order_err_ = 0;
    } else if (method_ == "rk45") {
      // Pair of order 5(4), first same as last, derived by a standalone tool
      a_ = {0, 0, 0, 0, 0, 0, 0,
            2./9, 0, 0, 0, 0, 0, 0,
            95./1024, 225./1024, 0, 0, 0, 0, 0,
            595749./800000, -558657./160000, 44649./12500, 0, 0, 0, 0,
            8133487./3056130, -547./42, 22178368./1898505, -4601000./12530133, 0, 0, 0,
            418551./146000, -404253./29200, 197614472./16087375, -74305./320251,
            -4991679./67174600, 0, 0,
            629./6600, 0, 677888./1454475, 208000./337799, -453789./920200, 73./231, 0};
      alpha_ = {0, 2./9, 5./16, 33./40, 20./21, 1, 1};
      m_ = {629./6600, 0, 677888./1454475, 208000./337799, -453789./920200, 73./231, 0};
      // Dense output of order 4, continuously differentiable: coefficients of theta, ..., theta^4
      n_dense_ = 4;
      m_dense_ = {1, 0, 0, 0, 0, 0, 0,
                  -3063./1100, 0, 1948928./484825, -1404000./337799, 2398599./460100,
                  -292./77, 3./2,
                  2434./825, 0, -8982016./1454475, 520000./48257, -1426194./115025,
                  292./33, -4,
                  -471./440, 0, 254208./96965, -2028000./337799, 1231713./184040,
                  -365./77, 5./2};
      // Error estimate: difference with the embedded solution of order 4
      m_err_ = {89834./99443025, 0, -183980032./58439351025, 2874688000./122151834189,
                -92439186./1540529825, 468368./9281349, -8./681};
      order_err_ = 4;
    } else {
      casadi_error("'erk': unknown method '" + method_ + "', expected 'rk45' or 'rk4'");
    }
  }

  void Erk::set_erk_prob() {
    // First same as last: the last stage is evaluated at t0 + h and y1
    casadi_int s = alpha_.size();
    pe_.fsal = alpha_[s - 1] == 1;
    for (casadi_int j = 0; j < s; ++j) pe_.fsal = pe_.fsal && a_[(s - 1) * s + j] == m_[j];
  }

  void Erk::set_work(void* mem, const double**& arg, double**& res,
      casadi_int*& iw, double*& w) const {
    auto m = static_cast<ErkMemory*>(mem);
    RuntimeIntegrator::set_work(mem, arg, res, iw, w);
    m->e.prob = &pe_;
    m->e.nom = nullptr;
    if (replay_) {
      m->e_nom.prob = &pe_;
      m->e.nom = &m->e_nom;
    }
  }

  int Erk::kernel(RuntimeIntegratorMemory* mem) const {
    auto m = static_cast<ErkMemory*>(mem);
    return casadi_erk(&m->d, &m->e);
  }

  void Erk::codegen_method(CodeGenerator& g) const {
    g.add_auxiliary(CodeGenerator::AUX_ERK);
    g << "pe.fsal = " << pe_.fsal << ";\n";
    g << "e.prob = &pe;\n";
    if (replay_) {
      g << "pen.fsal = " << pe_.fsal << ";\n";
      g << "en.prob = &pen;\n";
      g << "e.nom = &en;\n";
    } else {
      g << "e.nom = 0;\n";
    }
  }

  Erk::Erk(DeserializingStream& s) : RuntimeIntegrator(s) {
    s.version("Erk", 1);
    s.unpack("Erk::method", method_);
    set_coefficients();
    set_ivp_prob();
    set_nominal_prob();
    set_erk_prob();
  }

  void Erk::serialize_body(SerializingStream &s) const {
    RuntimeIntegrator::serialize_body(s);
    s.version("Erk", 1);
    s.pack("Erk::method", method_);
  }

} // namespace casadi
