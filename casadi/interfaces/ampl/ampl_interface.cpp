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


#include "ampl_interface.hpp"
#include "casadi/core/casadi_misc.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOL_AMPL_EXPORT
  casadi_register_nlpsol_ampl(Nlpsol::Plugin* plugin) {
    plugin->creator = AmplInterface::creator;
    plugin->name = "ampl";
    plugin->doc = AmplInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &AmplInterface::options_;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_AMPL_EXPORT casadi_load_nlpsol_ampl() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_ampl);
  }

  AmplInterface::AmplInterface(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }


  AmplInterface::~AmplInterface() {
    clear_mem();
  }

  Options AmplInterface::options_
  = {{&Nlpsol::options_},
     {{"solver",
       {OT_STRING,
        "AMPL solver binary"}}
     }
  };

  void AmplInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    // Set default options
    solver_ = "ipopt";

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="solver") {
        solver_ = op.first;
      }
    }

    // Extract the expressions
    casadi_assert(oracle().is_a("SXFunction"),
                  "Only SX supported currently.");
    vector<SX> xp = oracle().sx_in();
    vector<SX> fg = oracle()(xp);

    // Get x, p, f and g
    SX x = xp.at(NL_X);
    SX p = xp.at(NL_P);
    SX f = fg.at(NL_F);
    SX g = fg.at(NL_G);
    casadi_assert(p.is_empty(), "'p' currently not supported");

    // Names of the variables, constraints
    vector<string> x_name, g_name;
    for (int i=0; i<nx_; ++i) x_name.push_back("x[" + str(i) + "]");
    for (int i=0; i<ng_; ++i) g_name.push_back("g[" + str(i) + "]");
    int max_x_name = x_name.back().size();
    int max_g_name = g_name.empty() ? 0 : g_name.back().size();

    // Calculate the Jacobian, gradient
    SX jac_g = SX::jacobian(g, x);
    SX grad_f = gradient(f, x);

cout << "g = " << g << endl;
    // Extract the shared subexpressions
    vector<SX> ex = {f, g, jac_g, grad_f}, v, vdef;
    shared(ex, v, vdef);
    f = ex[0];
    g = ex[1];
    jac_g = ex[2];
    grad_f = ex[3];

    // Print problem:
    uout() << "minimize f subject to lbx<=x<=ubx, lbg<=g<=ubg\n"
           << "with: f = " << g << "\n"
           << " g = " << jac_g << "\n"
           << " grad_f = " << grad_f << "\n"
           << " jac_g = " << grad_f << "\n"
           << v << " = " << vdef << "\n";

    // Create .nl file
    fstream nl("casadi.nl", fstream::out);
    // Header
    nl << "g3 1 1 0\n";
    // Type of constraints
    nl << nx_ << " " // number of variables
       << ng_ << " " // number of constraints
       << 1 << " " // number of objectives
       << 0 << " " // number of ranges
       << 0 << " " // ?
       << 0 << "\n"; // number of logical constraints

    // Nonlinearity - assume all nonlinear for now TODO: Detect
    nl << ng_ << " "  // nonlinear constraints
       << 1 << "\n"; // nonlinear objectives

    // Network constraints
    nl << 0 << " " // nonlinear
       << 0 << "\n"; // linear

    // Nonlinear variables
    nl << nx_ << " " // in constraints
       << nx_ << " " // in objectives
       << nx_ << "\n"; // in both

    // Linear network ..
    nl << 0 << " " // .. variables ..
       << 0 << " " // .. arith ..
       << 0 << " " // .. functions ..
       << 0 << "\n"; // .. flags

    // Discrete variables
    nl << 0 << " " // binary
       << 0 << " " // integer
       << 0 << " " // nonlinear
       << 0 << " " // ?
       << 0 << " " // ?
       << 0 << "\n"; // ?

    // Nonzeros in the Jacobian, gradients
    nl << jac_g.nnz() << " " // nnz in Jacobian
       << nx_ << "\n"; // nnz in gradients

    // Maximum name length
    nl << max_x_name << " " // constraints
       << max_g_name << "\n"; // variables



    // Close file
    nl.close();
  }

  int AmplInterface::init_mem(void* mem) const {
    if (Nlpsol::init_mem(mem)) return 1;
    auto m = static_cast<AmplInterfaceMemory*>(mem);

    return 0;
  }

  void AmplInterface::set_work(void* mem, const double**& arg, double**& res,
                                   int*& iw, double*& w) const {
    auto m = static_cast<AmplInterfaceMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

  }

  void AmplInterface::solve(void* mem) const {
    auto m = static_cast<AmplInterfaceMemory*>(mem);

    // Check the provided inputs
    check_inputs(mem);
}

} // namespace casadi
