//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2011 - 2014 The MINOTAUR Team.
//

/// \file qpOASESEngine.h
/// \author Christian Kirches, Interdisciplinary Center for Scientific
///         Computing (IWR), Heidelberg University, GERMANY
/// \date 09/09/2011
/// Declare the class qpOASESEngine.
///
/// qpOASES is open-source and can be obtained under the LGPL license
/// from http://www.qpoases.org
///


#ifndef MINOTAURQPOASESENGINE_H
#define MINOTAURQPOASESENGINE_H

#include <boost/shared_ptr.hpp>

#include "QPEngine.h"
#include "WarmStart.h"

#include <qpOASES.hpp>

namespace Minotaur {

  class   Timer;
  class   Environment;
  class   Problem;
  class   Solution;
  typedef boost::shared_ptr<Environment> EnvPtr;
  typedef boost::shared_ptr<Problem> ProblemPtr;
  typedef boost::shared_ptr<Solution> SolutionPtr;


  /**
   * The qpOASES perspective of a QP in sparse row compressed format
   */
  class qpOASESData {
  public:
    qpOASES::SQProblem       *qp;

    qpOASES::SymSparseMat *H;
    qpOASES::SparseMatrixRow    *A;

    qpOASES::sparse_int_t *Hir;
    qpOASES::sparse_int_t *Hjc;
    qpOASES::real_t       *Hval;

    qpOASES::sparse_int_t *Air;
    qpOASES::sparse_int_t *Ajc;
    qpOASES::real_t       *Aval;

    size_t  n;
    size_t  m;
    double *g;
    double *lb;
    double *ub;
    double *lbA;
    double *ubA;

    qpOASESData (int n_, int m_, int nnzh_, int nnza_);
    ~qpOASESData ();
  };


  class qpOASESPoint {
  public:
    double *x;
    double *y;
    double *lb;
    double *ub;
    double *lbA;
    double *ubA;
    qpOASES::Bounds bounds;
    qpOASES::Constraints constraints;
  };


  class qpOASESWarmStart : public WarmStart
  {
  public:
    /// Return true if warm start information is initialized, false
    /// otherwise.
    virtual bool hasInfo() { return 1; };

    /// Write to an output stream
    virtual void write(std::ostream &) const {};

  public:
    qpOASESPoint pt;
  };

  /**
   * qpOASES call and timing statistics
   */
  class qpOASESStats {
  public:
    UInt   calls;     /// Total number of calls to solve.
    UInt   strCalls;  /// Calls to solve while strong branching.
    double time;      /// Sum of time taken in all calls to solve.
    double strTime;   /// Time taken in strong branching alone.
    double cTime;     /// Time taken in copying data for strong-branching.
    UInt   iters;     /// Sum of number of iterations in all calls.
    UInt   strIters;  /// Number of iterations in strong branching alone.

    UInt n_total;
    UInt n_cold;
    UInt n_strbr;
    UInt n_dive;
    UInt n_warm;

    UInt piv_total;
    UInt piv_cold;
    UInt piv_strbr;
    UInt piv_dive;
    UInt piv_warm;

    double time_total;
    double time_cold;
    double time_strbr;
    double time_dive;
    double time_warm;		
  };


  /**
   * qpOASESEngine is used to solve QP problems using the parametric
   * primal-dual active set code qpOASES.
   * qpOASES finds a KT point for the bounded QP problem
   *
   *      minimize    f(x) = x'Hx/2 + g'x
   *
   *      subject to  lb  <=   x <= ub
   *                  lbA <= A.x <= ubA
   *
   * where x and g are n-vectors, H is a symmetric n*n matrix, and A is an
   * m*n matrix. If H is also positive semi-definite then the KT point is a
   * global solution, else a stationary point.
   */

  class qpOASESEngine : public QPEngine {

  public:
    friend class Problem;

    /// Default constructor.
    qpOASESEngine ();

    /// Constructor using given environment options.
    qpOASESEngine (EnvPtr env);

    /// Return an empty qpOASESEngine pointer.
    EnginePtr emptyCopy();

    /// Method to read the problem and initialize qpOASES.
    void load(ProblemPtr problem);

    /// Method to unload the current problem
    void clear();

    /// Solve the problem that was loaded and report the status.
    EngineStatus solve();

    /// Report the solution value from the last solve.
    double getSolutionValue();

    /// Report the solution.
    ConstSolutionPtr getSolution();

    /// Report the status of the last solve.
    EngineStatus getStatus();

    // Change bound on a constraint.
    void changeBound(ConstraintPtr cons, BoundType lu, double new_val);

    // Change bound on a variable.
    void changeBound(VariablePtr var, BoundType lu, double new_val);

    // Implement Engine::changeBound(VariablePtr, double, double).
    void changeBound(VariablePtr var, double new_lb, double new_ub);

    // Implement Engine::addConstraint() */
    void addConstraint(ConstraintPtr) {};

    // change objective.
    void changeObj(FunctionPtr f, double cb);

    // Convert 'min f' to 'min -f'.
    void negateObj();

    // Implement Engine::changeConstraint().
    void changeConstraint(ConstraintPtr con, LinearFunctionPtr lf,
                          const double & ub);

    // Implement Engine::getWarmStart(). // NULL for now.
    ConstWarmStartPtr getWarmStart();

    // Implement Engine::getWarmStartCopy(). // NULL for now.
    WarmStartPtr getWarmStartCopy();

    // Implement Engine::loadFromWarmStart().
    void loadFromWarmStart(WarmStartPtr );

    // Implement Engine::setIterationLimit().
    void setIterationLimit(int limit);

    // Implement Engine::resetIterationLimit().
    void resetIterationLimit();

    // Implement Engine::enableStrBrSetup()
    void enableStrBrSetup();

    // Implement Engine::disableStrBrSetup()
    void disableStrBrSetup();

    // Write statistics.
    void writeStats();

    // get name.
    std::string getName() const;

    /// Destroy.
    ~qpOASESEngine();

  private:

    /// Allocate the data structures for qpOASES.
    void load_();

    /**
     * Actually call bqpd to solve a QP using a specific mode. It is called
     * after all data has been set.
     */
    void solve_ (int mode, double *f);

    /// Copy primal and dual values of the solution from bqpd.
    void storeSol_ (double f);

    /// Fill sparsity pattern and values of the gradients.
    void setGradient_ ();

    /// Fill sparsity pattern and values of the Hessian.
    void setHessian_ ();

    /// Set the intial point for solving the QP.
    void setInitialPoint_ ();

    /// Copy variable bounds from the problem.
    void setVarBounds_ ();

    /// Copy constraint bounds from the problem.
    void setConsBounds_ ();

    /// Environment.
    EnvPtr env_;

    /// Problem that is loaded, if any.
    ProblemPtr problem_;

    /// Solution found by the engine.
    SolutionPtr sol_;

    /// If a constraint is modified, this is set to true.
    bool consModed_;
    bool hessModed_;

    /// Holds the qpOASES view of a QP
    qpOASESData *data_;

    /// If true, we should try to resolve in a different mode when error is reported.
    bool resolveError_;

    /// Number of iterations that can be performed during solve
    int iterLimit_;

    /// Maximum number of iterations that can be performed during solve
    int maxIterLimit_;

    /// Bounds are considered infinite if their value goes beyond this value.
    const double infty_;

    /// Bounds are considered equalities if they are closer than this value.
    const double bTol_;

    /// Statistics
    qpOASESStats stats_;

    /// True if currently doing strong-branching iterations. False otherwise.
    bool strBr_;
    bool ldWarm_;
    UInt n_changed_;
    qpOASESWarmStart ws;
    double objOff_;

  private:
    void savePoint (qpOASESPoint *);
    void restorePoint (qpOASESPoint *);
  };

  typedef boost::shared_ptr<qpOASESEngine> qpOASESEnginePtr;
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
