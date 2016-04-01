//
//     MINOTAUR -- It's only 1/2 bull  
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//
/**
 * \file LinFeasPump.h
 * \brief Feasibility pump using Linear relaxation
 * \author Jayash Koshal, Argonne National Laborator
 */

#ifndef LINFEASPUMP_H
#define LINFEASPUMP_H

#include "FeasibilityPump.h"

namespace Minotaur {
  class Engine;
  class LinearFunction;
  class LinearHandler;
  class Problem;
  class QGHandler;
  typedef boost::shared_ptr<Engine> EnginePtr;
  typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;
  typedef boost::shared_ptr<LinearHandler> LinHandlerPtr;
  typedef boost::shared_ptr<Problem> ProblemPtr;
  typedef boost::shared_ptr <QGHandler> QGHandlerPtr;

  /// statistics for Linear Feasibility Pump
  struct LinFeasStats {
    double bestObjValue; /// Objective value for best feasible sol
    UInt numLPs;         /// Number of NLPs solved in the heuristic
  };

  /**
   * \brief Linear Feasibility Pump for MINLPs.
   *
   * A Linear Feasibility Pump heuristic used to find solutions for Mixed
   * Integer NLPs by solving a linear relaxation of NLP using a LP
   * engine. An NLP is solved after every "n" iterations and solution
   * from NLP solve is used to construct a LP relaxation. This class is
   * derived from the class FeasibilityPump. 
   */

  class LinFeasPump : public FeasibilityPump {

  public:

    /**
     * \brief Default constructor
     *
     * \param[in] env Environment pointer
     * \param[in] p Problem pointer
     * \param[in] e1 The NLP engine pointer
     * \param[in[ e2 The LP engine pointer
     *
     * The constructor initializes the base class using env, p and e1
     * and uses e2 to initialize the LP engine. Logger of the base is
     * also defined in this function.
     */
    LinFeasPump(EnvPtr env, ProblemPtr p, EnginePtr e1, EnginePtr e2);

    /// Default destructor
    ~LinFeasPump();

    /// Call to the heuristic
    void solve(NodePtr node, RelaxationPtr rel, SolutionPoolPtr s_pool);

    /// Write statistics to the logger
    void writeStats(std::ostream &out) const;

  protected:

    /// Message name for the heuristic
    const static std::string me_;

    /// gradient of the objective function
    double* gradientObj_;

    /// Linear Handler pointer
    LinHandlerPtr lh_;

    /// LP Engine to be used to solving linear relaxation
    EnginePtr lpE_;

    /// objective improvement constraint pointer
    ConstraintPtr objConstraint_;

    /// The objective variable added by linearization of objective. If
    /// the objective is linear, it is NULL.
    VariablePtr objVar_;

    /// clone of linear objective function
    LinearFunctionPtr olfClone_;

    /// QG Handler
    QGHandlerPtr qh_;

    /// Relaxation Pointer
    RelaxationPtr r_;

    /// Statistics
    LinFeasStats* statsLFP_;

    // base class method.
    void constructObj_(ProblemPtr prob, ConstSolutionPtr sol);

    /**
     * \brief Fucntion to implement the linear feasibility pump
     *
     * \param[in] s_pool Pointer to solution pool
     * \param[in] x Is not required. Can be NULL.
     */
    void implementFP_(const double* x, SolutionPoolPtr s_pool);

    /**
     * \brief Calculate the gap between the NLP relaxation solution and
     * integer feasible solution.
     *
     * \param[in] f_nlp nlp solution value.
     * \param[in] x const pointer to integer feasible solution.
     *
     * \return gap. gap is given by 
     * \f$ \frac{f(x)-f(x_nlp)}{|f(x)|+1e-6}\f$
     */
    double getSolGap_(double f_nlp, double f_feas);

    /** 
     * \brief A function to prepare the linear relaxation
     *
     * \return True if problem is infeasible else false
     *
     * This function first initializes the Linear Handler and QG
     * Handler which are then used to construct a linear relaxation
     * from the original problem
     */
    bool prepareLP_();

    /**
     * \brief This function makes a cut by including the objective as
     * constraint.
     *
     * \param[in] f_nlp NLP solution value
     * solution
     * \param[in] x const pointer to the primal feasible solution
     *
     * If a feasible solution is found we add an additional constraint
     * to the problem to force the heuristic to find a better feasible
     * fucntion. The constraint added is a linear relaxation of the
     * objective at the feasible point, i.e.
     * \f$ \nabla f(x^{LP}) + (x-x^{LP}) < 0 \f$
     */
    void separatingCut_(double f_nlp, SolutionPoolPtr s_pool);

    /**
     * \brief Function to decide whether to use Linear Feasibility Pump
     *
     * \return true if to be implemented else fasle
     *
     * The linear feasibility pump is not implemented for problems with
     * equality constraints for the form \f$ c(x) = 0 \f$
     */
    bool shouldFP_();

  };

  typedef boost::shared_ptr<LinFeasPump> LinFeasPumpPtr;
}
#endif

// local variables:
// mode: c++
// eval: (c-set-style "gnu")
// eval: (setq indent-tabs-mode nil)
// end:
