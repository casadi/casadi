// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

/**
 * \file FilterSQPEngine.h
 * \author Sven Leyffer, Argonne National Laboratory.
 * Define the class FilterSQPEngine.
 */


#ifndef MINOTAURFILTERSQPENGINE_H
#define MINOTAURFILTERSQPENGINE_H

#include "NLPEngine.h"
#include "WarmStart.h"

namespace Minotaur {
  class   Environment;
  class   FilterSQPWarmStart;
  class   Problem;
  class   Solution;
  class   Timer;
  typedef boost::shared_ptr<Environment> EnvPtr;
  typedef boost::shared_ptr<FilterSQPWarmStart> FilterWSPtr;
  typedef boost::shared_ptr<const FilterSQPWarmStart> ConstFilterWSPtr;
  typedef boost::shared_ptr<Problem> ProblemPtr;
  typedef boost::shared_ptr<Solution> SolutionPtr;

  struct FilterSQPStats {
    UInt calls;     ///<  Total number of calls to solve.
    UInt iters;     ///<  Sum of number of iterations in all calls. 
    UInt strCalls;  ///<  Calls to solve while strong branching.
    UInt strIters;  ///<  Number of iterations in strong branching alone.
    double strTime; ///<  time taken in strong branching alone.
    double time;    ///<  Sum of time taken in all calls to solve.
  };

  class FilterSQPWarmStart : public WarmStart {
  public:

    /// Default constructor
    FilterSQPWarmStart(); 

    /// Copy constructor. Creates a full copy, not just copies pointers.
    FilterSQPWarmStart(ConstFilterWSPtr warm_st);

    /// Destroy
    ~FilterSQPWarmStart();

    /// Get solution.
    SolutionPtr getPoint();

    // Implement WarmStart::hasInfo().
    bool hasInfo();

    /**
     * Overwrite the primal and dual values of warm-start. Sometimes, the
     * warm-start data is initialized and needs to be updated. This
     * should be called in place of deleting and creating a new warm-start
     * object.
     */
    void setPoint(SolutionPtr sol);

    // Implement WarmStart::write().
    void write(std::ostream &out) const;

  private:
    /// The starting solution that is used to warm-start.
    SolutionPtr sol_;
  };

  /**
   * FilterSQPEngine is the class that is used to solve NLP problems using
   * the FilterSQP solver. FilterSQP can be used to solve problems completely to
   * optimality or approximately.
   */
  class FilterSQPEngine : public NLPEngine {

  public:
    friend class Problem;

    /// Default constructor.
    FilterSQPEngine();    

    /// Constructor using given environment options.
    FilterSQPEngine(EnvPtr env);    

    /// Destroy.
    ~FilterSQPEngine();

    // Base class method
    void addConstraint(ConstraintPtr c);

    // Change bound on a constraint.
    void changeBound(ConstraintPtr cons, BoundType lu, double new_val);

    // Change bound on a variable.
    void changeBound(VariablePtr var, BoundType lu, double new_val);

    // Implement Engine::changeBound(VariablePtr, double, double).
    void changeBound(VariablePtr var, double new_lb, double new_ub);

    // Implement Engine::changeConstraint().
    void changeConstraint(ConstraintPtr con, LinearFunctionPtr lf, 
                          double lb, double ub);

    // Implement Engine::changeConstraint().
    void changeConstraint(ConstraintPtr c, NonlinearFunctionPtr nlf);

    // change objective.
    void changeObj(FunctionPtr f, double cb);

    void clear();

    /// Restore settings after strong branching.
    void disableStrBrSetup();

    /// Return an empty FilterSQPEngine pointer.
    EnginePtr emptyCopy();

    /// Make settings for strong branching.
    void enableStrBrSetup();

    /**
     * Evaluate the activity of constraints and fill the values in 'c'.
     * This is a callback function that is called by filter. 
     */
    void evalCons(const double *x, double *c, int *error);

    /**
     * Evaluate the gradient of objective and fill in 'a'. 
     * This is a callback function that is called by filter. 
     */
    void evalGrad(const double *x, double *a, int *error);

    /**
     * Evaluate the hessian of lagrangian at point x and fill in. 
     * This is a callback function that is called by filter. 
     */
    void evalHessian(const double *x, double *lam, 
                     const int phase, double *ws, int *lws, int *l_hess,
                     int *li_hess, int *error);

    /// Report the solution value from the last solve.
    double evalObjValue(const double *x, int *err);

    // Get name.
    std::string getName() const;

    /// Report the solution value from the last solve. 
    double getSolutionValue();

    /// Report the solution.
    ConstSolutionPtr getSolution();

    /// Report the status of the last solve.
    EngineStatus getStatus();

    // Implement Engine::getWarmStart(). // NULL for now.
    ConstWarmStartPtr getWarmStart();

    // Implement Engine::getWarmStartCopy(). // NULL for now.
    WarmStartPtr getWarmStartCopy();

    /// Method to read the problem and initialize FilterSQP.
    void load(ProblemPtr problem);

    // Implement Engine::loadFromWarmStart().
    void loadFromWarmStart(WarmStartPtr ws);

    // Convert 'min f' to 'min -f'.
    void negateObj();

    // base class method.
    void removeCons(std::vector<ConstraintPtr> &delcons);

    // Implement Engine::resetIterationLimit().
    void resetIterationLimit();

    // Implement Engine::setIterationLimit().
    void setIterationLimit(int limit);

    /// Solve the problem that was loaded and report the status.
    EngineStatus solve();

    // Write statistics.
    void writeStats(std::ostream &out) const;

  private:

    /// Jacobian storage.
    double *a_;

    /// Lower bounds.
    double *bl_;

    /// if lb > ub and lb < ub+bTol_, we make ub = lb.
    double bTol_;

    /// Upper bounds.
    double *bu_;

    /// values of constraint functions
    double *c_;

    /**
     * If true, reallocate space in the next solve. Important, if
     * constraints or objectives have changed.
     */
    bool consChanged_;

    /// Linear (L) or Nonlinear (N).
    char *cstype_;

    /// Environment.
    EnvPtr env_;

    /// if rstat_[4]<feasTol_, then the solution is considered feasible.
    const double feasTol_;

    /// Statistics.
    int *istat_;

    /// Number of iterations that can be performed during solve
    int iterLimit_; 

    /**
     * la_ stores the sparsity pattern of jacobian. It needs to be
     * evaluated only once.
     */
    int *la_;

    /// Lagrange multipliers.
    double *lam_;

    /**
     * lws_ stores the sparsity pattern of hessian of lagrangian. It needs
     * to be evaluated only once. In our implementation:
     * lws_[0]    = number of entries in the hessian of lagrangian.
     * We set phl = 0 (phl is a common block in filter-sqp),
     * phr        = 1,
     * lws_[phr]  = lws_[1] = row index of the first entry in the
     *              hessian,
     * lws_[phr+i]= lws_[1+i] = row index of the i-th entry in the
     *              hessian,
     * phc        = hess_nz + 1,
     * lws_[phc]  = lws_[1+hess_nz] = column index of the first entry 
     *              in the hessian,
     * lws_[phc+i]= lws_[1+hess_nz+i] = column index of the i-th entry 
     *              in the hessian.
     */
    int *lws_;

    /// Copy of lws_;
    int *lws2_;

    /// Max value of iterLimit_, when solving a relaxation
    const int maxIterLimit_;

    /// String name used in log messages.
    static const std::string me_;

    /// Need to multiply lagrange multipliers by -1 in callback. Storage.
    double *mlam_;

    /**
     * True if we want to save warm start information of the current
     * solution for the next solve. False, if no information needs to be
     * saved. Saving information does not mean that it will be used.
     * useWs_ flag must be on to use the warm-start information.
     */
    bool prepareWs_;

    /// Problem that is loaded, if any.
    ProblemPtr problem_;

    /// Statistics.
    double *rstat_; 

    /// scale factors.
    double *s_;

    /**
     * If true, then copy the solution from the last solve. Otherwise,
     * don't copy it.
     */
    bool saveSol_;

    /// Solution found by the engine. 
    SolutionPtr sol_;

    /// Statistics.
    FilterSQPStats *stats_;

    /// If true, we are currently in strong-branching mode. False otherwise.
    bool strBr_;

    /// Timer
    Timer *timer_;

    /**
     * True if we want to use warm-start information, either from a
     * previous solve or from a user provided structure.
     */
    bool useWs_;

    /// warm start information.
    FilterWSPtr warmSt_;

    double *ws_;

    /// Solution.
    double *x_;

    /// Free arrays used by filter-sqp.
    void freeStorage_();

    /// Copy bounds from problem into filter-sqp's arrays.
    void setBounds_();

    /// Allocate storage space for filter-sqp.
    void setStorage_(int mxwk, int maxa);

    /// Allocate space and fill sparsity pattern of the Jacobian and Hessian.
    void setStructure_();

  };

  typedef boost::shared_ptr<FilterSQPEngine> FilterSQPEnginePtr;
} // end namespace Minotaur 

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
