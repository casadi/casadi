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


#include "dummy.hpp"

#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/calculus.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOL_DUMMY_EXPORT
      casadi_register_nlpsol_dummy(Nlpsol::Plugin* plugin) {
    plugin->creator = Dummy::creator;
    plugin->name = "dummy";
    plugin->doc = Dummy::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Dummy::options_;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_DUMMY_EXPORT casadi_load_nlpsol_dummy() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_dummy);
  }

  Dummy::Dummy(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }

  Dummy::~Dummy() {
    uout() << "Solver class destroyed" << std::endl;
  }

  Options Dummy::options_
  = {{&Nlpsol::options_},
     {{"my_string",
       {OT_STRING,
        "Example for a string-typed option."}},
      {"my_int",
       {OT_INT,
        "Example for an int-typed option."}},
      {"my_double",
       {OT_DOUBLE,
        "Example for a double-typed option."}},
     }
  };

  void Dummy::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    uout() << "Solver class init" << std::endl;

    // Avoid the superclass fixing lam_x
    bound_consistency_ = false;

    // Default options
    my_string_ = "foo";
    my_int_ = 0;
    my_double_ = 1;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="my_string") {
        my_string_ = op.second.to_string();
      } else if (op.first=="my_int") {
        my_int_ = op.second;
      } else if (op.first=="my_double") {
        my_double_ = op.second;
      }
    }

    // List the callbacks that you need

    // Get/generate required functions
    create_function("nlp_fg", {"x", "p"}, {"f", "g"});
    // First order derivative information
    Function nlp_jac_fg = create_function("nlp_jac_fg", {"x", "p"},
                      {"f", "grad:f:x", "g", "jac:g:x"});

    // Get the sparsity of the Jacobian
    J_ = nlp_jac_fg.sparsity_out(3);

    // Inspect the sparsity: check the user guide for an explaination on sparsity
    uout() << "Sparsity of dg/dx:" << std::endl;
    J_.spy();

    // Allocate space in the work vector for the (dense) objective gradient
    alloc_w(nx_, true);

    // Allocate space in work vector for nonzeros of jacobian
    alloc_w(J_.nnz(), true);

  }

  void Dummy::set_work(void* mem, const double**& arg, double**& res,
                                casadi_int*& iw, double*& w) const {
    auto m = static_cast<DummyMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Set gf pointer to the first available location in the work vector
    m->gf = w; w += nx_;

    // Set Jk pointer to the first available location in the work vector
    m->Jk = w; w += J_.nnz();

  }

  int Dummy::solve(void* mem) const {
    auto m = static_cast<DummyMemory*>(mem);

    uout() << "Solve"  << std::endl;

    uout() << "Initial guess for 'x':" <<
      std::vector<double>(m->x, m->x+nx_) << std::endl;
    uout() << "Initial guess for 'lam_x':" <<
      std::vector<double>(m->lam_x, m->lam_x+nx_) << std::endl;
    uout() << "Initial guess for 'lam_g':" <<
      std::vector<double>(m->lam_g, m->lam_g+ng_) << std::endl;
    uout() << "Parameter values:" <<
      std::vector<double>(m->p, m->p+np_) << std::endl;


    // Set x,p inputs
    m->arg[0] = m->x;
    m->arg[1] = m->p;

    // Set f, grad:f:x, g, jac:g:x outputs
    m->res[0] = &m->f;
    m->res[1] = m->gf;
    m->res[2] = m->g;
    m->res[3] = m->Jk;

    // Compute
    calc_function(m, "nlp_jac_fg");

    uout() << "Objective gradient (always dense): " <<
      std::vector<double>(m->gf, m->gf+nx_) << std::endl;
    uout() << "Constraint jacobian (nonzeros, column-major): " <<
      std::vector<double>(m->Jk, m->Jk+J_.nnz())<< std::endl;

    // This is not an actual solver.
    // Let's just output some nonsense answer
    for (int i=0;i<nx_;++i)
      m->x[i]= 666+2*i;

    std::fill(m->lam_g, m->lam_g+ng_, 2.5);
    std::fill(m->lam_x, m->lam_x+nx_, 3.5);

    // Invent some statistics
    m->iter_count = 42;

    return 0;
  }

  Dict Dummy::get_stats(void* mem) const {
    Dict stats = Nlpsol::get_stats(mem);
    auto m = static_cast<DummyMemory*>(mem);
    stats["iter_count"] = m->iter_count;
    return stats;
  }
} // namespace casadi
