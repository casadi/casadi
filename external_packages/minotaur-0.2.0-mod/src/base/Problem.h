//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//


/**
 * \file Problem.h
 * \brief Declare base class Problem.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURPROBLEM_H
#define MINOTAURPROBLEM_H

#include<ios>
#include "Types.h"
#include "Variable.h"

namespace Minotaur {

  class Engine;
  class Function;
  class HessianOfLag;
  class Jacobian;
  class LinearFunction;
  class NonlinearFunction;
  class Objective;
  struct ProblemSize;
  class QuadraticFunction;
  class SOS;
  class SparseMatrix;
  typedef boost::shared_ptr<Function> FunctionPtr;
  typedef boost::shared_ptr<Jacobian> JacobianPtr;
  typedef boost::shared_ptr<HessianOfLag> HessianOfLagPtr;
  typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;
  typedef boost::shared_ptr<NonlinearFunction> NonlinearFunctionPtr;
  typedef boost::shared_ptr<Objective> ObjectivePtr;
  typedef boost::shared_ptr<ProblemSize> ProblemSizePtr;
  typedef boost::shared_ptr<QuadraticFunction> QuadraticFunctionPtr;
  typedef boost::shared_ptr<const ProblemSize> ConstProblemSizePtr;
  typedef SOS* SOSPtr;

  /**
   * \brief The Problem that needs to be solved.
   *
   * The Problem class contains the Variables, Constraints, Objectives that
   * collectively define a particular problem that can be solved. A
   * problem can be described as
   * 
   * min/max f(x)
   * 
   * s.t.  l_i <= g_i(x) <= u_i, i = 1,2, ... , m
   * 
   *       l_i <= x_i <= u_i, i = 1, 2, ..., n
   * 
   * A way of setting up an problem is to first ask it to create 'n'
   * variables. Then the constraints and objective can be added one by one.  
   * 
   * A Problem is a very generic class. The Relaxation classes
   * are derived from it. 
   */
  class Problem {
  public:
    /// Default constructor
    Problem();

    /// Destroy
    virtual ~Problem();

    /// Add 'c' to both lb and ub of a constraint.
    virtual void addToCons(ConstraintPtr cons, double c);

    /// Add a linear function to the objective.
    virtual void addToObj(LinearFunctionPtr lf);

    /// Add a constant term to the objective.
    virtual void addToObj(double cb);

    /// Fill up the statistics about the size of the problem into size_.
    virtual void calculateSize(bool shouldRedo=false);

    /// Change a bound (lower or upper) on a variable with ID=id.
    virtual void changeBound(UInt id, BoundType lu, double new_val);

    /// Change both bounds (lower and upper) on a variable with ID=id
    virtual void changeBound(UInt id, double new_lb, double new_ub);

    /// Change a bound (lower or upper) on a variable 'var'. 
    virtual void changeBound(VariablePtr var, BoundType lu, double new_val);

    /// Change lower and upper bounds on the variable 'var'
    virtual void changeBound(VariablePtr var, double new_lb, 
                             double new_ub); 

    /// Change a bound (lower or upper) on a constraint 'con'. 
    virtual void changeBound(ConstraintPtr con, BoundType lu, double new_val);

    /// Change lower and upper bounds on the constraint 'con'
    virtual void changeBound(ConstraintPtr con, double new_lb, 
                             double new_ub); 

    /**
     * \brief Change the linear function, and the bounds of a constraint.
     * \param [in] con Original constraint that is to be changed.
     * \param [in] lf The new linear function.
     * \param [in] lb The new lower bound.
     * \param [in] ub The new upper bound.
     */
    virtual void changeConstraint(ConstraintPtr con, LinearFunctionPtr lf, 
                                  double lb, double ub);

    /**
     * \brief Change the nonlinear function and bounds of a constraint.
     * \param [in] con Original constraint that is to be changed.
     * \param [in] nlf The new nonlinear function.
     */
    virtual void changeConstraint(ConstraintPtr con, NonlinearFunctionPtr nlf);


    /**
     * \brief Replace the objective function with a new function. 
     *
     * \param[in] f The new obejctive function. f is cloned. If f is modified
     * after this call, it won't affect the objective.
     * \param[in] cb The new objective constant.
     */
    virtual void changeObj(FunctionPtr f, double cb);

    /**
     * \brief Check whether variables used in the constraints belong to the
     * problem or not.
     *
     * This is a sanity check, used only for debugging.
     * \returns 1 if the check failed, 0 if passed.
     */
    virtual int checkConVars() const;

    /**
     * \brief Delete the whole Problem.
     *
     * Variables and constraints are so interlinked that we just can not call
     * the destructor. This function just deletes all the constraints. The
     * variables and functions can still be used after this is called.
     */
    virtual void clear();

    /**
     * \brief Clone the given Problem class. Jacobian and Hessian in the cloned
     * problem are NULL.
     *
     * The variables are created. If the functions are stored in native format,
     * they are also cloned. Problem size and the initial point are cloned as
     * well.
     */
    ProblemPtr clone() const;

    /// Iterate over constraints. Returns the 'begin' iterator.
    virtual ConstraintConstIterator consBegin() const 
    { return cons_.begin(); }

    /// Iterate over constraints. Returns the 'end' iterator.
    virtual ConstraintConstIterator consEnd() const { return cons_.end(); }

    /// Delete marked constraints.
    virtual void delMarkedCons();

    /// Delete marked variables.
    virtual void delMarkedVars();

    /**
     * \brief Return what type of problem it is. May result in re-calculation of
     * the problem size.
     */
    virtual ProblemType findType();

    /// Return a pointer to the constraint with a given index
    virtual ConstraintPtr getConstraint(UInt index) const;

    /// Return the hessian of the lagrangean. Could be NULL.
    virtual HessianOfLagPtr getHessian() const;

    /**
     * \brief Get the initial point. Used by some engines like IpoptEngine. The
     * pointer returned from this function should not be changed or deleted.
     */
    virtual const double * getInitialPoint() const { return initialPt_; }

    /// Return the jacobian. Could be NULL.
    virtual JacobianPtr getJacobian() const;

    /// Get pointer to the log manager. Could be NULL.
    virtual LoggerPtr getLogger();

    /// Return the number of constraints.
    virtual UInt getNumCons() const { return cons_.size(); }

    /// Return the number of constraints marked for deletion.
    virtual UInt getNumDCons() const { return numDCons_; }

    /// Return the number of variables marked for deletion.
    virtual UInt getNumDVars() const { return numDVars_; }

    /**
     * \brief Return the number of non-zeros in the hessian of the lagrangean of the 
     * problem.
     *
     * The lagrangean is defined as:
     * \\sigma . f(x) + \\sum_{i=0}^{m-1}\\lambda_i . g_i(x), 
     * where \\sigma \\in R^1 and \\lambda \\in R^m are the dual multipliers. 
     * The hessian, w.r.t. x, is thus a square symmetric matrix. usually the
     * multipliers are provided by NLP solvers. 
     * Such solvers may require during initialization, the number of non-zeros
     * in the lower triangular of the hessian.
     */
    virtual UInt getNumHessNnzs() const;

    /// Return the number of non zerors in the jacobian of the constraints.
    virtual UInt getNumJacNnzs() const;

    /// Return the number of linear constraints in the problem.
    UInt getNumLinCons(); 

    /// Return the number of SOS Type 1 constraints.
    UInt getNumSOS1();

    /// Return the number of SOS Type 2 constraints.
    UInt getNumSOS2();

    /// Return the number of variables.
    virtual UInt getNumVars() const { return vars_.size(); }

    /// Return a pointer to the objective Function
    virtual ObjectivePtr getObjective() const;

    /// Return the value of objective function at given point x.
    double getObjValue(const double *x, int *err) const;

    /// Fill up the statistics about the size of the problem into size_.
    ConstProblemSizePtr getSize() const;

    /// Return a pointer to the variable with a given index
    virtual VariablePtr getVariable(UInt index) const;

    /**
     * \brief Return true if the derivative is available through Minotaur's own
     * routines for storing nonlinear functions.
     */
    virtual bool hasNativeDer() const;

    /**
     * \brief Returns true if the problem has only linear constraints and linear
     * objectives.
     */
    virtual bool isLinear();

    /// Return true if a constraint is marked deleted.
    virtual bool isMarkedDel(ConstConstraintPtr con);

    /// Return true if a constraint is marked deleted.
    virtual bool isMarkedDel(ConstVariablePtr var);

    /**
     * \brief Returns true if the problem has 
     * (1) linear or quadratic objective, and
     * (2) linear constraints only.
     */
    virtual bool isQP();

    /**
     * \brief Returns true if the problem has only linear or quadratic constraints
     * and linear or quadratic objectives.  Returns false if a problem is
     * linear. Returns false if problem is nonlinear.
     */
    virtual bool isQuadratic();

    /**
     * \brief Mark a constraint for deleting.
     * 
     * The constraint is not deleted, just marked. Call Problem::delMarkedCons()
     * to actually delete all the marked constraints.
     * \param[in] con The constraint to be marked.
     */
    virtual void markDelete(ConstraintPtr con);

    /**
     * \brief Mark a variable as deleted.
     *
     * The variable is not deleted, just marked. Call Problem::delMarkedVars()
     * to actually delete all the marked variables.
     * \param[in] var The variable to be marked.
     */
    virtual void markDelete(VariablePtr var);

    /// The objective is multiplied by -1.
    virtual void negateObj();

    /**
     * \brief Add a new binary variable and return a pointer to it. A name is
     * automatically generated by default.
     */
    virtual VariablePtr newBinaryVariable(); 

    /**
     * \brief Add a new binary variable.
     *
     * \param[in] name The predefined name for this variable.
     */
    virtual VariablePtr newBinaryVariable(std::string name); 

    /**
     * \brief Add a new constraint and return a pointer to it. A name is
     * automatically generated by default.
     *
     * \param[in] f Pointer to the Function in the constraint. It is not cloned.
     * The pointer is saved as it is.
     * \param[in] lb The lower bound of the constraint. May be -INFINITY.
     * \param[in] ub The upper bound of the constraint. May be +INFINITY.
     */
    virtual ConstraintPtr newConstraint(FunctionPtr f, double lb, double ub);

    /**
     * \brief Add a new constraint and return a pointer to it.
     *
     * \param[in] f Pointer to the Function in the constraint. It is not cloned.
     * The pointer is saved as it is. 
     * \param[in] lb The lower bound of the constraint. May be -INFINITY.
     * \param[in] ub The upper bound of the constraint. May be +INFINITY.
     * \param[in] name The name for the constraint.
     */
    virtual ConstraintPtr newConstraint(FunctionPtr f, double lb, double ub, 
                                        std::string name);

    /**
     * \brief Add a new objective. A name is automatically generated by
     * default.
     *
     * \param[in] f Pointer to the Function in the objective. It is not
     * cloned. The pointer is saved as it is.
     * \param[in] cb The constant term in the objective function.
     * \param[in] otyp Whether the objective is to Minimize or Maximize.
     */
    virtual ObjectivePtr newObjective(FunctionPtr f, double cb, 
                                      ObjectiveType otyp);

    /** 
     * \brief Add a new objective. 
     *
     * \param[in] f Pointer to the Function in the objective. It is not cloned.
     * The pointer is saved as it is.
     * \param[in] cb The constant term in the objective function.
     * \param[in] otyp Whether the objective is to Minimize or Maximize.
     * \param[in] name The name for the objective function.
     *
     * \returns Pointer to the newly added objective function.
     */
    virtual ObjectivePtr newObjective(FunctionPtr f, double cb, 
                                      ObjectiveType otyp, std::string name);

    /** 
     * \brief Add a new SOS constraint with a name. 
     *
     * \param[in] n Number of variables in this SOS constraint.
     * \param[in] type SOS1 (SOS type 1) or SOS2 (SOS type 2).
     * \param[in] weights Values of coefficients of variables in the SOS
     * constraint or just relative weights.
     * \param[in] vars Variables in the constraint.
     * \param[in] priority The priority provided by the user for this
     * constraint.
     * \param[in] name The name provided by the user for this
     * SOS.
     *
     * \returns Pointer to the newly added SOS data.
     */
    virtual SOSPtr newSOS(int n, SOSType type, const double *weights,
                          const VarVector &vars, int priority, std::string name);

    /** 
     * \brief Add a new SOS constraint (name generated automatically). 
     *
     * \param[in] n Number of variables in this SOS constraint.
     * \param[in] type SOS1 (SOS type 1) or SOS2 (SOS type 2).
     * \param[in] weights Values of coefficients of variables in the SOS
     * constraint or just relative weights.
     * \param[in] vars Variables in the constraint.
     * \param[in] priority The priority provided by the user for this
     * constraint.
     *
     * \returns Pointer to the newly added SOS data.
     */
    virtual SOSPtr newSOS(int n, SOSType type, const double *weights,
                          const VarVector &vars, int priority);


    /**
     * \brief Add a new continuous, unbounded variable to the Problem. 
     * \param[in] stype The source of the variable
     */
    virtual VariablePtr newVariable(VarSrcType stype=VarOrig); 

    /**
     * \brief Add a new variable with bounds, type. A name is automatically
     * generated by default.
     *
     * \param[in] lb The lower bound on the variable. May be -INFINITY.
     * \param[in] ub The upper bound on the variable. May be +INFINITY.
     * \param[in] vtype Type of the variable: Integer, Continuous, Binary.
     * \param[in] stype The source of the variable
     */
    virtual VariablePtr newVariable(double lb, double ub, VariableType vtype,
                                    VarSrcType=VarOrig); 

    /**
     * \brief Add a new variable.
     *
     * \param[in] lb The lower bound on the variable. May be -INFINITY.
     * \param[in] ub The upper bound on the variable. May be +INFINITY.
     * \param[in] vtype Type of the variable: Integer, Continuous, Binary.
     * \param[in] name Name of the variable.
     * \param[in] stype The source of the variable
     */
    virtual VariablePtr newVariable(double lb, double ub, VariableType vtype,
                                    std::string name, VarSrcType=VarOrig); 

    /**
     * \brief Clone the variables pointed by the iterators and add them.
     *
     * Given starting and stopping iterators of variables, clone these
     * variables and add the copies to this problem. Do not add them to any
     * constraints or objectives. The IDs are not copied.
     * \param[in] v_begin The 'begin' iterator of the variable vector.
     * \param[in] v_end The 'end' iterator of the variable vector.
     * \param[in] stype The source of the variables
     */
    virtual void newVariables(VariableConstIterator v_begin, 
                              VariableConstIterator v_end,
                              VarSrcType stype=VarOrig);

    /**
     * \brief Setup problem data-structures for solving it.
     *
     * Prepare necessary data structures for the next solve. e.g.
     * if constraints have been modified, then re-evalate the sparsity
     * pattern of Jacobian and Hessian.
     */
    virtual void prepareForSolve();

    /// Remove objective from the Problem.
    virtual void removeObjective(); 

    /// Remove the quadratic part of objective and return it.
    virtual QuadraticFunctionPtr removeQuadFromObj();

    /**
     * Remove the jacobian and hessian data structures. Useful when you want to
     * re-compute the derivatives after a problem has been modified.
     */
    virtual void resetDer();

    /**
     * \brief Reverse the sense of a constraint.
     * 
     * \param[in] cons The constraint whose sense has to be reversed.
     */
    virtual void reverseSense(ConstraintPtr cons);

    /**
     * \brief Set the engine that is used to solve this problem.
     *
     * The problem contains a pointer to the engine so that whenever the problem
     * is modified, the engine also gets the modifications. This function sets
     * the engine that must be modified whenever the problem is modified.
     * \param[in] engine The engine pointer.
     */
    virtual void setEngine(Engine* engine);

    /**
     * \brief Set an initial point.
     *
     * Initial point is used by some engines like IpoptEngine. If the
     * initial point has already been set before, it is overwritten by the
     * new point.
     * \param[in] x An array of double values containing the coordinates of the
     * initial point.
     */
    virtual void setInitialPoint(const double *x);

    /** 
     * \brief Set an initial point.
     *
     * Same as function Problem::setInitialPoint, but only set values for the
     * first 'k' variables. Put in because of AMPL's defined variables.
     * \param[in] x An array of double values containing the coordinates of the
     * initial point.
     * \param[in] k The first 'k' variables will be initialized.
     */
    virtual void setInitialPoint(const double *x, size_t k);

    /**
     * \brief Add a pointer to the hessian of the Lagrangean. 
     *
     * \param[in] hessian Pointer to the HessianOfLag object.
     */
    virtual void setHessian(HessianOfLagPtr hessian);

    /**
     * \brief Set the jacobian of the constraints. 
     *
     * \param[in] jacobian Pointer to the Jacobian object.
     */
    virtual void setJacobian(JacobianPtr jacobian);

    /**
     * Set the log manager.
     *
     * \param[in] logger The log manager which should be used for logging.
     */
    virtual void  setLogger(LoggerPtr logger);

    /**
     * \brief Ask Problem to construct its own jacobian and hessian using
     * Minotaur's native code for nonlinear functions.
     */
    void setNativeDer();

    /**
     * \brief Change the variable type.
     *
     * \param[in] var The variable pointer whose type needs to be changed.
     * \param[in] type The new VariableType.
     */
    virtual void setVarType(VariablePtr var, VariableType type);

    virtual SOSConstIterator sos1Begin() const { return sos1_.begin(); };
    virtual SOSConstIterator sos1End() const { return sos1_.end(); };
    virtual SOSConstIterator sos2Begin() const { return sos2_.begin(); };
    virtual SOSConstIterator sos2End() const { return sos2_.end(); };

    /**
     * \brief Substitute a variable 'out' with the variable 'in' through out the
     * problem.
     *
     * \param[in] out The variable that is to be substituted out.
     * \param[in] in The variable that replaces the variable 'out'.
     * \param[in] rat The ratio of substitution.
     * \f$v_{in} = rat \times v_{out}\f$.
     */
    virtual void subst(VariablePtr out, VariablePtr in, double rat=1.0);

    /// Should be called in the Engine's destructor
    virtual void unsetEngine();

    /// Iterate over variables.
    virtual VariableConstIterator varsBegin() const { return vars_.begin(); }

    /// Iterate over variables.
    virtual VariableConstIterator varsEnd() const { return vars_.end(); }

    /// only for debugging, developing etc.
    virtual void write(std::ostream &out, std::streamsize out_p=6) const;

    /// Write the problem size to logger_
    virtual void writeSize(std::ostream &out) const;

  protected:
    /// Vector of constraints.
    ConstraintVector cons_;

    /**
     * \brief Flag that is turned on if the constraints are added or modified.
     *
     * When the problem is changed, say constraints are added or deleted,
     * all associated changes do not happen at that time. For example, when
     * new constraint is added, the hessian and jacobian do not change.
     * This variable is true if constraints have changed since the last
     * time all changes were applied.
     */
    bool consModed_;

    /// Engine that must be updated if problem is loaded to it, could be null 
    Engine* engine_;

    /// Pointer to the hessian of the lagrangean. Could be NULL.
    HessianOfLagPtr hessian_;

    /// Initial point. Can be NULL.
    double * initialPt_;

    /// Pointer to the jacobian of constraints. Can be NULL.
    JacobianPtr jacobian_;

    /// Pointer to the log manager. All output messages are sent to it.
    LoggerPtr logger_;

    /// For logging
    static const std::string me_;

    /// If true, set up our own Hessian and Jacobian.
    bool nativeDer_;

    /// ID of the next constraint.
    UInt nextCId_;

    /// ID of the next SOS constraint.
    int nextSId_;

    /// ID of the next variable.
    UInt nextVId_;

    /// Number of constraints marked for deletion
    UInt numDCons_;

    /// Number of variables marked for deletion
    UInt numDVars_;	

    /// Objective, could be NULL.
    ObjectivePtr obj_;

    /// Size statistics for this Problem.
    ProblemSizePtr size_;

    /// SOS1 constraints.
    SOSVector sos1_;

    /// SOS2 constraints.
    SOSVector sos2_;

    /// Vector of variables.
    VarVector vars_;

    /// True if variables delete, added or their bounds changed.
    bool varsModed_;

    /// Count the types of constraints and fill the values in size_.
    virtual void countConsTypes_();

    /// Count the types of objectives and fill the values in size_.
    virtual void countObjTypes_();

    /// Fill up the size_ with number of variables of each type.
    virtual void countVarTypes_();

    /**
     * \brief Update the function types of all variables based on whether they
     * appear as linear, quadratic or nonlinear in each objective and
     * constraint.
     */
    virtual void findVarFunTypes_();

    bool isPolyp_();

    void setIndex_(VariablePtr v, UInt i);

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
