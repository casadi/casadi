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

#ifndef IPOPT_INTERNAL_HPP
#define IPOPT_INTERNAL_HPP

#include "ipopt_solver.hpp"
#include "casadi/fx/nlp_solver_internal.hpp"

/** \brief  Forward declarations */
namespace Ipopt{
  struct IpoptApplication;
}

namespace CasADi{

/**
@copydoc NLPSolver_doc
*/
class IpoptInternal : public NLPSolverInternal{
friend class IpoptUserClass;

public:
  explicit IpoptInternal(const FX& F, const FX& G, const FX& H, const FX& J, const FX& GF);
  virtual ~IpoptInternal();
  virtual IpoptInternal* clone() const{ return new IpoptInternal(*this);}
  
  FX getGF() const { return GF_; }
  
virtual void init();
virtual void evaluate(int nfdir, int nadir);

protected:
  
Ipopt::IpoptApplication* app;
void *userclass;
std::map<std::string,opt_type> ops_;

  // The NLP functions
  /// Gradient of the objective function
  FX GF_; 

  // Ipopt callback functions
  bool eval_f(int n, const double* x, bool new_x, double& obj_value);
  bool eval_grad_f(int n, const double* x, bool new_x, double* grad_f);
  bool eval_g(int n, const double* x, bool new_x, int m, double* g);
  bool eval_jac_g(int n, const double* x, bool new_x,int m, int nele_jac, int* iRow, int *jCol,double* values);
  bool eval_h(const double* x, bool new_x, double obj_factor, const double* lambda,bool new_lambda, int nele_hess, int* iRow,int* jCol, double* values);
  void finalize_solution(const double* x, const double* z_L, const double* z_U, const double* g, const double* lambda, double obj_value);
  bool get_bounds_info(int n, double* x_l, double* x_u,int m, double* g_l, double* g_u);
  bool get_starting_point(int n, bool init_x, double* x,bool init_z, double* z_L, double* z_U,int m, bool init_lambda,double* lambda);
  void get_nlp_info(int& n, int& m, int& nnz_jac_g,int& nnz_h_lag);
  int get_number_of_nonlinear_variables() const;
  bool get_list_of_nonlinear_variables(int num_nonlin_vars, int* pos_nonlin_vars) const;
  bool intermediate_callback(const double* x, const double* z_L, const double* z_U, const double* g, const double* lambda, double obj_value, int iter, double inf_pr, double inf_du,double mu,double d_norm,double regularization_size,double alpha_du,double alpha_pr,int ls_trials);
  
  // Accummulated time since last reset:
  double t_eval_f_; // time spent in eval_f
  double t_eval_grad_f_; // time spent in eval_grad_f
  double t_eval_g_; // time spent in eval_g
  double t_eval_jac_g_; // time spent in eval_jac_g
  double t_eval_h_; // time spent in eval_h
  double t_callback_fun_;  // time spent in callback function
  double t_callback_prepare_; // time spent in callback preparation
  
};

} // namespace CasADi

#endif //IPOPT_INTERNAL_HPP
