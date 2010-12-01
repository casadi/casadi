#ifndef IPOPT_INTERNAL_HPP
#define IPOPT_INTERNAL_HPP

#include "ipopt_solver.hpp"
#include "casadi/fx/nlp_solver_internal.hpp"

/** \brief  Forward declarations */
namespace Ipopt{
  struct IpoptApplication;
}

namespace CasADi{
    
class IpoptInternal : public NLPSolverInternal{
friend class IpoptUserClass;

public:
  explicit IpoptInternal(const FX& F, const FX& G, const FX& H, const FX& J);
  virtual ~IpoptInternal();

virtual void init();
virtual void evaluate(int fsens_order, int asens_order);

protected:
  
Ipopt::IpoptApplication* app;
void *userclass;
std::map<std::string,opt_type> ops_;



/** \brief  options */
  bool verbose_;           // verbosity
  bool exact_hessian_; // use exact hessian

  bool eval_f(int n, const double* x, bool new_x, double& obj_value);
  bool eval_grad_f(int n, const double* x, bool new_x, double* grad_f);
  bool eval_g(int n, const double* x, bool new_x, int m, double* g);
  bool eval_jac_g(int n, const double* x, bool new_x,int m, int nele_jac, int* iRow, int *jCol,double* values);
  bool eval_h(const double* x, bool new_x, double obj_factor, const double* lambda,bool new_lambda, int nele_hess, int* iRow,int* jCol, double* values);

  void finalize_solution(const double* x, const double* z_L, const double* z_U, const double* g, const double* lambda, double obj_value);
  bool get_bounds_info(int n, double* x_l, double* x_u,int m, double* g_l, double* g_u);
  bool get_starting_point(int n, bool init_x, double* x,bool init_z, double* z_L, double* z_U,int m, bool init_lambda,double* lambda);
  void get_nlp_info(int& n, int& m, int& nnz_jac_g,int& nnz_h_lag);


};

} // namespace CasADi

#endif //IPOPT_INTERNAL_HPP
