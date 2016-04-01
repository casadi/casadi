// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

/// \file BqpdEngine.h
/// \author Sven Leyffer, Argonne National Laboratory.
/// Declare the class BqpdEngine.


#ifndef MINOTAURBQPDENGINE_H
#define MINOTAURBQPDENGINE_H

#include <boost/shared_ptr.hpp>

#include "QPEngine.h"

namespace Minotaur {

  class   Timer;
  class   Environment;
  class   Problem;
  class   Solution;
  typedef boost::shared_ptr<Environment> EnvPtr;
  typedef boost::shared_ptr<Problem> ProblemPtr;
  typedef boost::shared_ptr<Solution> SolutionPtr;

  struct BqpdStats {
    UInt calls;     /// Total number of calls to solve.
    UInt strCalls;  /// Calls to solve while strong branching.
    double time;    /// Sum of time taken in all calls to solve.
    double strTime; /// Time taken in strong branching alone.
    double cTime;   /// Time taken in copying data for strong-branching.
    UInt iters;     /// Sum of number of iterations in all calls. 
    UInt strIters;  /// Number of iterations in strong branching alone.
  };

  class BqpdData;
  /**
   * BqpdEngine is used to solve QP problems using bqpd. 
   * bqpd finds a KT point for the bounded QP problem
   * 
   *      minimize    f(x) = ct.x + xt.G.x/2
   * 
   *      subject to  l <= [I : A]t.x <= u                  (t = transpose)
   * 
   * where x and c are n-vectors, G is a symmetric n*n matrix, and A is an
   * n*m matrix. If G is also positive semi-definite then the KT point is a
   * global solution, else usually a local solution. The method may also be
   * used efficiently to solve an LP problem (G=0). bqpd can be used to solve
   * problems completely to optimality or approximately.
   */

  class BqpdEngine : public QPEngine {
  public:
    friend class Problem;

    /// Default constructor.
    BqpdEngine();    

    /// Constructor using given environment options.
    BqpdEngine(EnvPtr env);    

    /// Destroy.
    ~BqpdEngine();

    // Implement Engine::addConstraint() */
    void addConstraint(ConstraintPtr);

    // Change bound on a constraint.
    void changeBound(ConstraintPtr cons, BoundType lu, double new_val);

    // Change bound on a variable.
    void changeBound(VariablePtr var, BoundType lu, double new_val);

    // Implement Engine::changeBound(VariablePtr, double, double).
    void changeBound(VariablePtr var, double new_lb, double new_ub);

    // Implement Engine::changeConstraint().
    void changeConstraint(ConstraintPtr con, LinearFunctionPtr lf, 
                          double lb, double ub);

    // base class method
    void changeConstraint(ConstraintPtr con, NonlinearFunctionPtr nlf);

    // change objective.
    void changeObj(FunctionPtr f, double cb);

    /// Method to unload the current problem
    void clear();

    // Implement Engine::disableStrBrSetup()
    void disableStrBrSetup();

    /// Return an empty BqpdEngine pointer.
    EnginePtr emptyCopy();

    // Implement Engine::enableStrBrSetup()
    void enableStrBrSetup();

    // get name.
    std::string getName() const;

    /// Report the solution.
    ConstSolutionPtr getSolution();

    /// Report the solution value from the last solve. 
    double getSolutionValue();

    /// Report the status of the last solve.
    EngineStatus getStatus();

    // Implement Engine::getWarmStart(). // NULL for now.
    ConstWarmStartPtr getWarmStart() {return WarmStartPtr(); }; 

    // Implement Engine::getWarmStartCopy(). // NULL for now.
    WarmStartPtr getWarmStartCopy() {return WarmStartPtr(); };

    /// Method to read the problem and initialize bqpd.
    void load(ProblemPtr problem);

    // Implement Engine::loadFromWarmStart().
    void loadFromWarmStart(WarmStartPtr ) {};

    // Convert 'min f' to 'min -f'.
    void negateObj();

    // delete constraints.
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
    /// Number of bound changes since last solve.
    UInt bndChanges_;

    /**
     * \brief True if some variable bounds are relaxed. Bqpd has difficulty
     * hot-strating in this case.
     */
    bool bndRelaxed_;

    const double bTol_;

    /// Checkpoint copy of fStart_.
    BqpdData *chkPt_;

    /// If a constraint is modified, this is set to true. 
    bool consModed_;

    /// Array to calculate dual solution for constraints.
    double *dualCons_;

    /// Array to calculate dual solution for variables.
    double *dualX_;

    /// Environment.
    EnvPtr env_;

    /// Information for full start
    BqpdData *fStart_;

    /// Bounds are considered infinite if their value goes beyond this value.
    const double infty_;

    /// Number of iterations that can be performed during solve
    int iterLimit_; 

    /// Max value of iterLimit_, when solving a relaxation
    const int maxIterLimit_;

    /// String name used in log messages.
    static const std::string me_;

    /// Constant part of the obj.
    double objOff_;

    /// If the previous call was a strong-branching call.
    bool prevStrBr_;

    /// Problem that is loaded, if any.
    ProblemPtr problem_;

    /**
     * If true, we should try to resolve in a different mode when error is
     * reported.
     */
    bool resolveError_;

    /// Solution found by the engine. 
    SolutionPtr sol_;

    /// Statistics.
    BqpdStats *stats_;

    /// True if currently doing strong-branching iterations. False otherwise.
    bool strBr_;

    /// Timer for bqpd solves.
    Timer *timer_;

    /// Mode used for warm starting: 1-6
    int wsMode_;

    /// Free the memory allocated
    void freeProb_();

    /// Allocate the data structures for Bqpd.
    void load_();

    /// Copy constraint bounds from the problem.
    void setConsBounds_();

    /// Fill sparsity pattern and values of the gradients.
    void setGradient_();

    /// Fill sparsity pattern and values of the Hessian.
    void setHessian_();

    /// Set the intial point for solving the QP.
    void setInitialPoint_();

    /// Copy variable bounds from the problem.
    void setVarBounds_();

    /**
     * Actually call bqpd to solve a QP using a specific mode. It is called
     * after all data has been set.
     */
    void solve_(int mode, double &f);

    /// Copy primal and dual values of the solution from bqpd.
    void storeSol_(double f);

  };


  /// Information for restarting from the previous optimal solution.
  class BqpdData {
  public: 
    /// Constructor.
    BqpdData(UInt n_t, UInt m_t, int kmax_t, UInt maxa_t, UInt lh1_t,
             UInt nJac, bool zero=true);

    /// Destroy.
    ~BqpdData();

    /// Allocate space and copy.
    BqpdData *clone();

    /// Only copy. No space allocation.
    void copyFrom(const BqpdData* rhs);

    /// Display all data.
    void write(std::ostream &out) const;

    /// Number of variables.
    UInt n;

    /// Number of constraints.
    UInt m;

    /// kmax given to bqpd
    int kmax;

    /// Number of nonzeros in Hessian.
    UInt lh1;

    /// Number of nonzeros in Jacobian.
    UInt nJac;

    /// Size of a
    UInt maxa;

    /// Initial point for solving QP.
    double *x;

    /// Residuals/multipliers.
    double *r;

    /// Steepest-edge normalization coefficients .
    double *e;

    /// Denominators for ratio tests.
    double *w;

    /// Gradient vector of f(x).
    double *g;

    /// Indices of the active constraints .
    int *ls;

    /// Workspace associated with recursion.
    double *alp;

    /// Workspace associated with recursion.
    int *lp;

    /// Information on return from bqpd.
    int *info;

    /// Lower bounds for variables and constraints.
    double *bl;

    /// Upper bounds for variables and constraints.
    double *bu;

    /// Storage for jacobian.
    double *a;

    /// Storage for jacobian.
    int *la;

    /// Storage for hessian values and other things.
    double *ws;

    /// Storage for hessian indices and other things.
    int *lws;

    /// Pointer to equality constraints, used by bqpd.
    int peq;

    // Dimension of reduced-space, set only when mode>=0.
    int k;
  };

  
  typedef boost::shared_ptr<BqpdEngine> BqpdEnginePtr;
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
