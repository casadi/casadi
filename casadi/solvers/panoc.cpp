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


#include "panoc.hpp"

#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/calculus.hpp"

// include the C panoc solver
#include "panoc.h"

#ifndef SUCCESS
  #define SUCCESS 0
#endif // !SUCCESS

#ifndef FAILURE
  #define FAILURE 1
#endif // !FAILURE

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOL_PANOC_EXPORT
      casadi_register_nlpsol_panoc(Nlpsol::Plugin* plugin) {
    plugin->creator = Panoc::creator;
    plugin->name = "panoc";
    plugin->doc = Panoc::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Panoc::options_;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_PANOC_EXPORT casadi_load_nlpsol_panoc() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_panoc);
  }

  Panoc::Panoc(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }

  Panoc::~Panoc() {
    uout() << "Solver class destroyed" << std::endl;
  }

/* types : OT_INT,OT_DOUBLE,OT_STRING */
  Options Panoc::options_
  = {{&Nlpsol::options_},
     {{"outer_iterations",
       {OT_INT,
        "Maximum number of loops solving the lagrangian."}},
      {"inter_interations",
       {OT_INT,
        "Maximum number of loops solving with PANOC."}},
      {"tolerance",
       {OT_DOUBLE,
        "Tolerance of PANOC solver"}},
      {"constraints_weight",
       {OT_DOUBLE,
        "weight of the contraints"}},
     }
  };

  void Panoc::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    uout() << "Solver class init" << std::endl;

    // Avoid the superclass fixing lam_x
    bound_consistency_ = false;

    // Default options
    this->outer_iterations=10;
    this->inner_iterations=300;
    this->tolerance=0.001;
    this->constraints_weight=100.;

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

    // First order derivative information
    Function nlp_jac_fg = create_function("nlp_jac_fg", {"x", "p"},
                      {"f", "grad:f:x");


    create_function("nlp_g", {"x", "p"}, {"g"});
    create_function("nlp_fwd_g", {"x", "p","fwd:x"},{"fwd:g"});

    // Allocate space in the work vector for the (dense) objective gradient
    alloc_w(nx_, true);

    // Allocate space in work vector for nonzeros of jacobian
    alloc_w(J_.nnz(), true);

  }

  void Panoc::set_work(void* mem, const double**& arg, double**& res,
                                casadi_int*& iw, double*& w) const {
    auto m = static_cast<DummyMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Set gf pointer to the first available location in the work vector
    m->gf = w; w += nx_;

    // Set Jk pointer to the first available location in the work vector
    m->Jk = w; w += J_.nnz();

  }

  int Panoc::solve(void* mem) const {
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

    // // Set x,p inputs
    // m->arg[0] = m->x;
    // m->arg[1] = m->p;

    // // Set f, grad:f:x, g, jac:g:x outputs
    // m->res[0] = &m->f;
    // m->res[1] = m->gf;
    // m->res[2] = m->g;
    // m->res[3] = m->Jk;

    struct optimizer_extended_problem extended_problem;

    this->outer_iterations=10;
    this->inner_iterations=300;
    this->tolerance=0.001;
    this->constraints_weight=100.;

    extended_problem.max_loops=this->outer_iterations;
    extended_problem.problem.dimension=2;
    extended_problem.problem.solver_params.tolerance=this->tolerance;
    extended_problem.problem.solver_params.max_interations=this->inner_iterations;
    extended_problem.problem.solver_params.buffer_size=15;

    extended_problem.constraints = this.constraints;
    extended_problem.constraints_forwad_diff = this.constraints_forward_diff;

    double lower_bound=-1;
    double upper_bound=1;
    double constraint_weight=1;

    double solution[2]={0,0};

    optimizer_init_extended_box(&extended_problem,lower_bound,upper_bound,constraint_weight);
    solve_extended_problem(&solution);
    optimizer_cleanup_extended();

    // Compute
    // calc_function(m, "nlp_jac_fg");

    // uout() << "Objective gradient (always dense): " <<
    //   std::vector<double>(m->gf, m->gf+nx_) << std::endl;
    // uout() << "Constraint jacobian (nonzeros, column-major): " <<
    //   std::vector<double>(m->Jk, m->Jk+J_.nnz())<< std::endl;

    // // This is not an actual solver.
    // // Let's just output some nonsense answer
    // for (int i=0;i<nx_;++i)
    //   m->x[i]= 666+2*i;

    // std::fill(m->lam_g, m->lam_g+ng_, 2.5);
    // std::fill(m->lam_x, m->lam_x+nx_, 3.5);

    // // Invent some statistics
    // m->iter_count = 42;

    return 0;
  }

  Dict Panoc::get_stats(void* mem) const {
    Dict stats = Nlpsol::get_stats(mem);
    auto m = static_cast<DummyMemory*>(mem);
    stats["iter_count"] = m->iter_count;
    return stats;
  }
} // namespace casadi

int Panoc::constraints(const double* x,double* out){
  out[0] = cos(x[0]+x[1])+0.5;
  out[1] = sin(x[0]) + 0.5;

  return SUCCESS;
}

int Panoc::constraints_forward_diff(const double* x,const double* y,double* out){
  out[0]=-sin(x[0]+x[1])*y[0]+ -sin(x[0]+x[1])*y[1];
  out[1]=cos(x[0])*y[0]+0*y[1];

  return SUCCESS;
}

double Panoc::rosen_cost_gradient(const double* input,double* gradient){
  gradient[0]=input[1]*2;
  gradient[1]= (2*sinh(input[1]))/(cosh(input[1])*cosh(input[1])*cosh(input[1]));

  return input[0]*input[0] + tanh(input[1])*tanh(input[1]);
}