// 
//   MINOTAUR -- It's only 1/2 bull
// 
//   (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file Engine.h
 * \brief Define the base class Engine.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURENGINE_H
#define MINOTAURENGINE_H

#include "Types.h"

namespace Minotaur {

  class   Constraint;
  class   Function;
  class   LinearFunction;
  class   NonlinearFunction;
  class   Solution;
  class   WarmStart;
  class   Engine;
  typedef boost::shared_ptr<Engine> EnginePtr;
  typedef boost::shared_ptr<const Engine> ConstEnginePtr;  
  typedef boost::shared_ptr<Function> FunctionPtr;
  typedef boost::shared_ptr<Constraint> ConstraintPtr;
  typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;
  typedef boost::shared_ptr<NonlinearFunction> NonlinearFunctionPtr;
  typedef boost::shared_ptr<const Solution> ConstSolutionPtr;
  typedef boost::shared_ptr<WarmStart> WarmStartPtr;
  typedef boost::shared_ptr<const WarmStart> ConstWarmStartPtr;

  /**
   * An Engine is a solver that can solve a Problem. In most invocations of an
   * engine, we will solve a Relaxation. This is an abstract base class and all
   * engines are implemented are derived from here.
   */
  class Engine {

  public:
    /// Default constructor.
    Engine();

    /// Destroy.
    virtual ~Engine();

    /// Add a new constraint to the engine.
    virtual void addConstraint(ConstraintPtr) = 0;

    /// Change a bound of a constraint. 
    virtual void changeBound(ConstraintPtr cons, BoundType lu, 
                             double new_val) = 0;

    /// Change a bound of a variable. 
    virtual void changeBound(VariablePtr var, BoundType lu, double new_val) 
      = 0;

    /// Change both bounds of a variable.
    virtual void changeBound(VariablePtr var, double new_lb, double new_ub) 
      = 0;

    /**
     * \brief Change the linear function, and the bounds of a constraint.
     * \param [in] c Original constraint that is to be changed.
     * \param [lf] The new linear function.
     * \param [lb] The new lower bound.
     * \param [ub] The new upper bound.
     */
    virtual void changeConstraint(ConstraintPtr c, LinearFunctionPtr lf, 
                                  double lb, double ub) = 0;

    /**
     * \brief Change the nonlinear function, and the bounds of a constraint.
     * \param [in] c Original constraint that is to be changed.
     * \param [nlf] The new nonlinear function.
     */
    virtual void changeConstraint(ConstraintPtr c, NonlinearFunctionPtr nlf) = 0;

    /// Change objective function.
    virtual void changeObj(FunctionPtr f, double cb) = 0;

    /// Clear the loaded problem, if any, from the engine.
    virtual void clear() = 0;

    /// Restore settings after strong branching.
    virtual void disableStrBrSetup() = 0;

    /// Get a fresh copy of the engine, without the problem loaded into it.
    virtual EnginePtr emptyCopy() {return EnginePtr();} // NULL by default.

    /// Make settings for strong branching.
    virtual void enableStrBrSetup() = 0;

    /// Get the solution obtained after solving the problem.
    virtual ConstSolutionPtr getSolution() = 0;

    /// Get the solution value.
    virtual double getSolutionValue() = 0;

    /// Solve the problem that was loaded previously.
    virtual EngineStatus solve() = 0;

    /// Get the name.
    virtual std::string getName() const = 0;

    /// Get the status of the last solve command.
    virtual EngineStatus getStatus() = 0;

    /// Return a string that describes the status in simple words.
    virtual std::string getStatusString();

    /**
     * Get warm start information from the engine. This warm start
     * information can change if the engine is used to solve something else
     * again.
     */
    virtual ConstWarmStartPtr getWarmStart() = 0;

    /**
     * Get a full copy of warm start information from the engine. Does not
     * change even if the engine starts solving something else later on.
     */
    virtual WarmStartPtr getWarmStartCopy() = 0;

    /**
     * Initialize the engine to solve the given problem. Memory is
     * allocated in this function.
     */
    virtual void load(ProblemPtr problem) = 0;

    /**
     * Use warm start information for solving the next problem. 
     * May Create a copy of WarmStart and use the
     * copy inside the engine; the copy (but not the original) gets updated
     * after solving a relaxation.
     */
    virtual void loadFromWarmStart(const WarmStartPtr ws) = 0;

    /// Negate the objective function. Min f is changed to Min -f.
    virtual void negateObj() = 0;

    /// Return pointer to the log manager
    virtual LoggerPtr  getLogger() { return logger_; }

    /**
     * \brief Delete constraints from the engine.
     *
     * \param [in] delcons A vector of constraint pointers that should be
     * deleted from the engine.
     */
    virtual void removeCons(std::vector<ConstraintPtr> &delcons) = 0;

    /// Reset the iteration limit to maximum possible.
    virtual void resetIterationLimit() = 0;

    /** 
     * Set a limit on number of iterations. For strong-branching, for
     * instance.
     */
    virtual void setIterationLimit(int limit) = 0;

    /// Set a new log manager
    virtual void setLogger(LoggerPtr logger) { logger_ = logger; }

    /// Set options to solve the NLP only once or very few times, with
    /// possibly several changes.
    virtual void setOptionsForSingleSolve() {};

    /// Set options to solve the NLP repeatedly, with few changes.
    virtual void setOptionsForRepeatedSolve() {};

    /**
     * Write statistics to the logger. If the log level is too low, no
     * statistics may be written.
     */
    virtual void writeStats(std::ostream &) const {};

  protected:
    /// Status of the last solve.
    EngineStatus status_;

    /// Keep log.
    LoggerPtr logger_;

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
