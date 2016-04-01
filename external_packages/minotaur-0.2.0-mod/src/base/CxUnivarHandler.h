//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file CxUnivarHandler.h
 * \brief Define the CxUnivarHandler class for handling convex univariate 
 * functions.
 * \author Ashutosh Mahajan, Argonne National Laboratory and Jim Luedtke,
 * UW-Madison
 */

#ifndef MINOTAURCXUNIVARHANDLER_H
#define MINOTAURCXUNIVARHANDLER_H

#include "Handler.h"

namespace Minotaur {

  class Engine;
  class Function;
  class LinearFunction;
  class Objective;
  class Problem;
  typedef boost::shared_ptr<Engine> EnginePtr;
  typedef boost::shared_ptr<Function> FunctionPtr;
  typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;
  typedef boost::shared_ptr<Objective> ObjectivePtr;
  typedef boost::shared_ptr<const Problem> ConstProblemPtr;


  class CxUnivarConstraintData {

  protected:

    /// Tolerance for constraint violation.
    double eTol_;

    /// Tolerance for when upper and lower bounds considered equal.
    double vTol_;

    /// This is the constraint of the orginal problem
    ConstraintPtr con_;

    /// Input variable, i.e., x in y = f(x)
    /// This should be a pointer to the `original' variable, look at
    /// corresponding bounds in the relaxation to get bounds
    ConstVariablePtr iv_;

    // Points to the relaxation version of the input variable
    // Updated only after relaxInit is called
    VariablePtr riv_;

    /// Output variable, i.e., y in y = f(x)
    /// This is a pointer to the `original' variable
    ConstVariablePtr ov_;

    // Points to the relaxation version of the output variable
    // Updated only after relaxInit is called
    VariablePtr rov_;

    /// 'L' => y <= f(x), so only add secants
    /// 'G' -> y >= f(x), so only add linearizations
    /// 'E' -> y == f(x), so do both
    char sense_;

    /// Secant constraint in the relaxation
    ConstraintPtr secCon_; 

    /// Array of linearization constraints in the relaxation
    ConstraintVector linCons_;

  public:

    /// Creates initial relaxations   
    void initRelax(RelaxationPtr rel, DoubleVector& tmpX, DoubleVector& grad);

    /// Update the current relaxation based on current variable bounds
    void updateRelax(RelaxationPtr rel, DoubleVector& tmpX, DoubleVector& grad,
                     ModVector &mods);

    bool isFeasible(const double* x);

    double getViol(const std::vector< double > & x);

    // Creates and adds the secant inequality defined by current constraint to
    // the relaxation
    // Returns the constraint that was added, which may be null if it was not
    // possible to add a constraint (e.g., due to an infinite bound)
    void addSecant(RelaxationPtr rel, ConstVariablePtr iv,
                   ConstVariablePtr ov, FunctionPtr fn, DoubleVector& tmpX,
                   bool init, ModVector &mods);

    // Creates and adds linearization inequalities to approximate the lower
    // envelope of the convex function
    // Returns a vector of constraints that were added
    void addLin(RelaxationPtr rel, ConstVariablePtr iv,
                ConstVariablePtr ov, FunctionPtr fn,
                DoubleVector& tmpX, DoubleVector& grad, bool init,
                ModVector &mods);

    ConstraintPtr getOriginalCon() const { return con_; }
    ConstraintPtr getSecantCon() const { return secCon_; }

    ConstVariablePtr getROutVar() const { return rov_; } ;

    ConstVariablePtr getRInputVar() const { return riv_; } ;

    //ConstraintVector & getLinCons
    ConstraintIterator linConsBegin()  { return linCons_.begin(); }
    ConstraintIterator linConsEnd() { return linCons_.end(); }


    char getSense() { return sense_; } ;

    /// Default constructor.
    CxUnivarConstraintData(double eTol, double vTol, ConstraintPtr newcon, ConstVariablePtr ovar,
                           ConstVariablePtr ivar, char sense);

    /// Destroy
    ~CxUnivarConstraintData() {};


  };

  typedef boost::shared_ptr<CxUnivarConstraintData> CxUnivarConstraintDataPtr;
  typedef std::vector<CxUnivarConstraintDataPtr> CxUnivarConstraintDataVector;
  typedef CxUnivarConstraintDataVector::iterator CxUnivarConstraintIterator;


  /**
   * An CxUnivarHandler handles convex univariate functions. The upper relaxation is
   * handled with the single secant inequality, the lower relaxation is handled
   * with linearizations.
   */

  /**
   *  TODO: This class could easily be extended to handle perspective cuts if a binary varaible is known which turns off the input variable
   *
   */
  class CxUnivarHandler : public Handler {


  protected:

    /// Tolerance for constraint violation.
    double eTol_;

    /// Tolerance for when upper and lower bounds considered equal.
    double vTol_;

    /// Original problem.
    ProblemPtr problem_;

    /// Logger.
    LoggerPtr logger_;

    /// For printing.
    static const std::string me_;

    /// Internal data associated with each constraint
    CxUnivarConstraintDataVector cons_data_;

    /// A temporary vector of zeros, for evaluating functions
    DoubleVector tmpX_;

    /// A temporary vector of zeros, for getting gradients
    DoubleVector grad_;

  public:
    /// Default constructor.
    CxUnivarHandler(EnvPtr env, ProblemPtr problem);

    /// Destroy
    ~CxUnivarHandler();

    /**
     * Adds constraint to list (as all handlers), but also constructs the
     * associated constraint data.
     */
    void addConstraint(ConstraintPtr newcon, ConstVariablePtr ivar,
                       ConstVariablePtr ovar, char sense = 'E'); 

    // base class method.
    void addConstraint(ConstraintPtr ) { assert(0); };

    /**
     *  For this handler, nothing is different at root or any node when doing full
     *  relax
     */
    void relaxInitFull(RelaxationPtr /*rel*/, bool* /* is_inf */) {};

    void relaxInitInc(RelaxationPtr rel, bool* is_inf);  

    /**
     * Check feasibility.
     */
    bool isFeasible(ConstSolutionPtr sol, RelaxationPtr relaxation, 
                    bool &should_prune, double &inf_meas);

    /**
     * Not implemented yet. Eventually, could add violated linearization
     * inequalities for underestimator portion
     */
    void separate(ConstSolutionPtr sol, NodePtr node, RelaxationPtr rel,
                  CutManager *cutman, SolutionPoolPtr s_pool, bool *sol_found,
                  SeparationStatus *status);


    /** 
     * Create a relaxation by adding the secant inequality for the upper estimator,
     * and some number of linearization inequalities, at a minimum from the end
     * points
     */
    virtual void relaxNodeFull(NodePtr /* node */, RelaxationPtr /* rel */,
                               bool* /* should_prune */) {} ; 

    /** 
     * Create a relaxation by updating the secant inequality for the upper
     * estimator, and adding lineariations at the end points, if they are new 
     */
    virtual void relaxNodeInc(NodePtr  node, RelaxationPtr rel,
                              bool* isInfeasible);

    // base class method
    virtual void getBranchingCandidates(RelaxationPtr rel, const DoubleVector &x,
                                        ModVector &mods, BrVarCandSet &cands,
                                        BrCandVector &gencands, bool & is_inf);

    // Implement Handler::getBrMod().
    virtual ModificationPtr getBrMod(BrCandPtr cand, DoubleVector &x, 
                                     RelaxationPtr rel, BranchDirection dir);

    // Implement Handler::getBranches().
    virtual Branches getBranches(BrCandPtr cand, DoubleVector & x,
                                 RelaxationPtr rel, SolutionPoolPtr s_pool);

    // presolve.
    virtual SolveStatus presolve(PreModQ *pre_mods, bool *changed);

    // Implement Handler::presolveNode().
    virtual bool presolveNode(RelaxationPtr p, NodePtr node,
                              SolutionPoolPtr s_pool, ModVector &p_mods,
                              ModVector &r_mods);


    // Write name
    virtual std::string getName() const;

  private:
    // Helper functions
    BranchPtr doBranch_(BranchDirection UpOrDown, ConstVariablePtr v,
                        double bvalue);



  };

  /// Shared pointer to CxUnivarHandler.
  typedef boost::shared_ptr<CxUnivarConstraintData> CxUnivarConstraintDataPtr;

  /// Shared pointer to const CxUnivarHandler.
  typedef boost::shared_ptr<const CxUnivarConstraintData>
    CxUnivarConstConstraintDataPtr;

  /// Shared pointer to CxUnivarHandler.
  typedef boost::shared_ptr<CxUnivarHandler> CxUnivarHandlerPtr;

  /// Shared pointer to const CxUnivarHandler.
  typedef boost::shared_ptr<const CxUnivarHandler> CxUnivarConstHandlerPtr;
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
