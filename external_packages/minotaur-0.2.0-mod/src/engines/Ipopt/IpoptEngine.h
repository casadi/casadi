// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file IpoptEngine.h
 * \brief Declare an interface to the IPOPT nlp-solver.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 * 
 * This file contains declartions for interface to Ipopt.
 */


#ifndef MINOTAURIPOPTENGINE_H
#define MINOTAURIPOPTENGINE_H

#include "NLPEngine.h"
#include "Solution.h"
#include "WarmStart.h"


namespace Ipopt {
  class IpoptApplication;
  class IpoptFunInterface;
  class TNLP;
}


namespace Minotaur {

  class   Environment;
  class   IpoptEngine;
  class   Problem;
  class   Timer;
  typedef boost::shared_ptr<Environment> EnvPtr;
  typedef boost::shared_ptr<IpoptEngine> IpoptEnginePtr;
  typedef boost::shared_ptr<Problem> ProblemPtr;

  struct IpoptStats {
    UInt calls;     ///< Total number of calls to solve.
    UInt iters;     ///< Sum of number of iterations in all calls. 
    double ptime;   ///< Sum of time taken in all calls to presolve.
    UInt strCalls;  ///< Calls to solve while strong branching.
    UInt strIters;  ///< Number of iterations in strong branching alone.
    double strTime; ///< time taken in strong branching alone.
    double time;    ///< Sum of time taken in all calls to solve.
  };


  /// The solution obtained from/for Ipopt.
  class   IpoptSolution;
  typedef boost::shared_ptr<IpoptSolution> IpoptSolPtr;
  typedef boost::shared_ptr<const IpoptSolution> ConstIpoptSolPtr;

  class IpoptSolution : public Solution {
  public:

    /// Default constructor
    IpoptSolution();

    /// Construct a solution for a problem from an array of doubles.
    IpoptSolution(const double *x, ProblemPtr problem); 

    /**
     * Construct a solution for a problem from an array of doubles and
     * objective value.
     */
    IpoptSolution(const double *x, double objval, ProblemPtr problem);

    /// Copy constructor
    IpoptSolution(ConstIpoptSolPtr);

    /// Destroy.
    ~IpoptSolution();

    /**
     * Return a pointer to the dual associated with the lower bound on
     * variables.
     */
    const double * getLowerDualOfVars() const {return dualXLow_;};

    /**
     * Return a pointer to the dual associated with the upper bound on
     * variables.
     */
    const double * getUpperDualOfVars() const {return dualXUp_;};

    // base class
    void setDualOfVars(const double *) { assert(!"implement me!"); };

    /// Set duals of variable-bounds: both the lower and upper.
    void setDualOfVars(const double *lower, const double *upper);

    /// Write to an output.
    void write(std::ostream &out) const;

  private:
    /// dual of lower bounds.
    double *dualXLow_;

    /// dual of upper bounds.
    double *dualXUp_;
  };


  class IpoptWarmStart;
  typedef boost::shared_ptr<IpoptWarmStart> IpoptWarmStartPtr;
  typedef boost::shared_ptr<const IpoptWarmStart> ConstIpoptWarmStartPtr;

  /// Class for saving and using Warm-start information in Ipopt.
  class IpoptWarmStart : public WarmStart {
  public:
    /// Default constructor
    IpoptWarmStart();

    /// Copy constructor. Creates a full copy, not just copies pointers.
    IpoptWarmStart(ConstIpoptWarmStartPtr ws);

    /// Destroy
    ~IpoptWarmStart();

    /// Return the soluton that can be used as starting point.
    IpoptSolPtr getPoint();

    // Implement WarmStart::hasInfo().
    bool hasInfo();

    void makeCopy();

    /**
     * Overwrite the primal and dual values of warm-start. Sometimes, the
     * warm-start data is initialized and needs to be updated. This
     * should be called in place of deleting and creating a new warm-start
     * object.
     */
    void setPoint(IpoptSolPtr sol);

    // Implement WarmStart::write().
    void write(std::ostream &out) const;

  private:
    /// The starting solution that is used to warm-start.
    IpoptSolPtr sol_;
  };


// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //

  class IpoptEngine : public NLPEngine {
  public:

    friend class Problem;

    /// Default constructor
    IpoptEngine();

    /// Default constructor
    IpoptEngine(EnvPtr env);

    /// Destroy
    ~IpoptEngine();

    // Base class method
    void addConstraint(ConstraintPtr c);

    void changeBound(ConstraintPtr cons, BoundType lu, double new_val);

    // Implement Engine::changeBound(VariablePtr, BoundType, double).
    void changeBound(VariablePtr var, BoundType lu, double new_val);

    // Implement Engine::changeBound(VariablePtr, double, double).
    void changeBound(VariablePtr var, double new_lb, double new_ub);

    // Implement Engine::changeConstraint().
    void changeConstraint(ConstraintPtr con, LinearFunctionPtr lf, 
                          double lb, double ub);

    // Implement Engine::changeConstraint().
    void changeConstraint(ConstraintPtr, NonlinearFunctionPtr);

    // change objective.
    void changeObj(FunctionPtr f, double cb);

    /// Clear the loaded problem, if any, from the engine.
    void clear();

    /// Restore settings after strong branching.
    void disableStrBrSetup();

    /// Return an empty IpoptEngine pointer.
    EnginePtr emptyCopy();

    /// Make settings for strong branching.
    void enableStrBrSetup();

    // get name.
    std::string getName() const;

    // Implement Engine::getSolution().
    ConstSolutionPtr getSolution();

    // Implement Engine::getSolutionValue().
    double getSolutionValue();

    // Implement Engine::getStatus().
    EngineStatus getStatus();

    // Implement Engine::getWarmStart().
    ConstWarmStartPtr getWarmStart();

    // Implement Engine::getWarmStartCopy().
    WarmStartPtr getWarmStartCopy();

    /// Load the problem into IPOPT. We create the TNLP interface to IPOPT.
    void load(ProblemPtr problem);

    // Implement Engine::loadFromWarmStart().
    void loadFromWarmStart(const WarmStartPtr ws);

    // Convert 'min f' to 'min -f'.
    void negateObj();

    // base class method.
    void removeCons(std::vector<ConstraintPtr> &delcons);

    // Implement Engine::resetIterationLimit().
    void resetIterationLimit();

    // Implement Engine::setIterationLimit().
    void setIterationLimit(int limit);

    void setOptionsForSingleSolve();

    void setOptionsForRepeatedSolve();

    /**
     * Solve the problem that was loaded. Uses the TNLP interface to IPOPT
     * for solving the NLP.
     */
    EngineStatus solve();

    /// Write statistics.
    void writeStats(std::ostream &out) const;

  private:
    /// True if a bound on variable has been modified since last change. 
    bool bndChanged_;

    /// True if a constraint has been modified since last change. 
    bool consChanged_;

    /// Optimal dual values associated with x >= l.
    double *dual_x_l_;

    /// Optimal dual values associated with x <= u.
    double *dual_x_u_;

    /// Optimal dual values associated with l <= g(x) <= u.
    double *dual_g_;

    /// Environment.
    EnvPtr env_;

    /// If lb>ub+etol_, declare infeasible.
    const double etol_;

    /// Where to put logs.
    LoggerPtr logger_;

    /// The maximum limit that can be set on Ipopt. 
    static const int maxIterLimit_ = 3000;

    /// String name used in log messages.
    static const std::string me_;

    Ipopt::IpoptApplication  *myapp_; 

    Ipopt::IpoptFunInterface *mynlp_;

    /**
     * True if we want to save warm start information of the current
     * solution for the next solve. False, if no information needs to be
     * saved. Saving information does not mean that it will be used.
     * useWs_ flag must be on to use the warm-start information.
     */
    bool prepareWs_;

    /// Try to solve without calling ipopt.
    bool presolve_();

    /// Problem that is loaded, if any.
    ProblemPtr problem_;

    /**
     * The solution obtained from last solve. The solution_ in 
     * IpoptFunInterface points to the same array as this one. So the
     * solution from Ipopt is copied to this array only once. This array
     * does not get freed in IpoptFunInterface and so must be freed here.
     */ 
    IpoptSolPtr sol_;

    /// Statistics.
    IpoptStats *stats_;

    /// True if strong branching. False otherwise.
    bool strBr_;

    /**
     * Timer for Ipopt solves. Includes time spent in call-backs:
     * evaluating functions and derivatives.
     */
    Timer *timer_;

    /**
     * True if we want to use warm-start information, either from a
     * previous solve or from a user provided structure.
     */
    bool useWs_;

    /// Private copy of the warm start information.
    IpoptWarmStartPtr ws_;

    /// Set problem specific options to make IPOPT faster
    void setOptionsForProb_();
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
