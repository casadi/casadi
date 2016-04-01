//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file IpoptEngineTnlp.h
 * \brief Declare the IPOPT TNLP class that provides Ipopt with the gradient,
 * hessian and function evaluations.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURIPOPTENGINETNLP_H
#define MINOTAURIPOPTENGINETNLP_H

namespace Ipopt {
  class IpoptFunInterface : public TNLP
  {
  public:

    /// Default constructor.
    IpoptFunInterface(Minotaur::ProblemPtr problem, 
                      Minotaur::IpoptSolPtr sol);

    /// default destructor.
    ~IpoptFunInterface();

    /// Method to return the objective value.
    bool eval_f(Index n, const Number* x, bool new_x, 
                Number& obj_value);

    /// Method to return the constraint residuals.
    bool eval_g(Index n, const Number* x, bool new_x, 
                Index m, Number* g);

    /// Method to return the gradient of the objective.
    bool eval_grad_f(Index n, const Number* x, 
                     bool new_x, Number* grad_f);

    /** 
     * Method to return:
     *   1) The structure of the jacobian (if "values" is NULL)
     *   2) The values of the jacobian (if "values" is not NULL).
     */
    bool eval_jac_g(Index n, const Number* x, 
                    bool new_x, Index m, Index nele_jac,
                    Index* iRow, Index *jCol,
                    Number* values);

    /** 
     * Method to return:
     *   1) The structure of the hessian of the lagrangian (if "values" is
     *   NULL) 
     *   2) The values of the hessian of the lagrangian (if "values"
     *   is not NULL)
     */
    bool eval_h(Index n, const Number* x, bool new_x,
                Number obj_factor, Index m, 
                const Number* lambda, bool new_lambda, 
                Index nele_hess, Index* iRow, 
                Index* jCol, Number* values);


    /** 
     * This method is called when the algorithm is complete so the TNLP can
     * store/write the solution 
     */
    void finalize_solution(SolverReturn status, 
                           Index n, const Number* x, const Number* z_L, 
                           const Number* z_U, Index m, const Number* g, 
                           const Number* lambda, Number obj_value, 
                           const IpoptData* ip_data, 
                           IpoptCalculatedQuantities* ip_cq);

    /// Method to return the bounds for my problem.
    bool get_bounds_info(Index n, Number* x_l,
                         Number* x_u, Index m, 
                         Number* g_l, Number* g_u);

    /// Method to return some info about the nlp.
    bool get_nlp_info(Index& n, Index& m, 
                      Index& nnz_jac_g, Index&
                      nnz_h_lag, IndexStyleEnum& index_style);

    /// Get solution.
    const Minotaur::IpoptSolPtr getSolution() {return sol_;}

    /// Get solution value.
    double getSolutionValue() const;

    /// Method to return the starting point for the algorithm.
    bool get_starting_point(Index n, bool init_x, 
                            Number* x, bool init_z, 
                            Number* z_L, Number* z_U,
                            Index m, bool init_lambda,
                            Number* lambda);

    /// Set solution.
    void setSolution(Minotaur::IpoptSolPtr sol) {sol_ = sol;}


  private:
    /// Copying is not allowed.
    IpoptFunInterface(const IpoptFunInterface &);

    /// Assignment is not allowed.
    IpoptFunInterface& operator=(const IpoptFunInterface&);

    /// If the variable is fixed, relax the upper bound by this much.
    double bOff_;

    /// If fabs(lb-ub)<bTol_ for a given variable, it is fixed to lb.
    double bTol_;

    /// Problem that is being solved.
    Minotaur::ProblemPtr problem_;

    /**
     * The final solution obtained from Ipopt. It is not clear from the
     * documentation of Ipopt that the final solution it sends in
     * finalize_solution above is not freed immediately after that function.
     * Hence we need to copy the solution.
     */
    Minotaur::IpoptSolPtr sol_;
  };
}
#endif
// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
