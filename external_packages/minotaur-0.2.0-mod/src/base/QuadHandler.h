//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file QuadHandler.h
 * \brief Declare a handler for simple quadratic constraints of the form
 * \f$ y_1 = x_1x_2 \f$,
 * and
 * \f$ y_1 = x_1^2 \f$,
 * It does not handle any other quadratic constraints.
 * \author Ashutosh Mahajan, IIT Bombay
 */

#ifndef MINOTAURQUADRATICHANDLER_H
#define MINOTAURQUADRATICHANDLER_H

#include "Handler.h"
#include "LinBil.h"

namespace Minotaur {

class LinearFunction;
class Timer;
typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;

/**
 * \brief A structure to save information about constraints of the form \f$ y
 * \leq x^2 \f$.
 */
struct LinSqr {
  VariablePtr   y;     ///> The variable y.
  VariablePtr   x;     ///> The variable x.
  ConstraintPtr oeCon; ///> The linear constraint that gives the over estimator
};
typedef LinSqr* LinSqrPtr;                 ///> Pointer to LinSqr
typedef std::vector<LinSqrPtr> LinSqrVec;  ///> Vector of LinSqr
typedef LinSqrVec::iterator LinSqrVecIter; ///> Iterator for LinSqr

/// Map of 'x' and the LinSqr that is used for \f$ y = x^2 \f$.
typedef std::map<VariablePtr, LinSqrPtr, CompareVariablePtr> LinSqrMap;
typedef LinSqrMap::iterator LinSqrMapIter; ///> Iterator for LinSqrMap


/**
 * A QuadHandler handles the quadratic functions of a problem in a simplistic
 * fashion. For now, we will just handle squares of singleton variables e.g.
 * \f$x_1^2 = y_1\f$, and bilinear terms: \f$x_1x_2 = y_1\f$.
 */
class QuadHandler : public Handler {
public:
  /// Default constructor.
  QuadHandler(EnvPtr env, ProblemPtr problem);

  /// Destroy
  ~QuadHandler();

  // base class method
  void addConstraint(ConstraintPtr newcon);

  // Implement Handler::getBranches().
  Branches getBranches(BrCandPtr cand, DoubleVector & x,
                       RelaxationPtr rel, SolutionPoolPtr s_pool);

  // base class method
  void getBranchingCandidates(RelaxationPtr rel, const DoubleVector &x,
                              ModVector &mods, BrVarCandSet &cands, 
                              BrCandVector &gencands, bool &is_inf);

  // base class method
  ModificationPtr getBrMod(BrCandPtr cand, DoubleVector &x, 
                           RelaxationPtr rel, BranchDirection dir);

  // base class method
  std::string getName() const;

  // base class method.
  bool isFeasible(ConstSolutionPtr sol, RelaxationPtr relaxation, 
                  bool &should_prune, double &inf_meas);

  // base class method.
  SolveStatus presolve(PreModQ *pre_mods, bool *changed);

  // base class method. Tightens bounds.
  bool presolveNode(RelaxationPtr p, NodePtr node,
                    SolutionPoolPtr s_pool, ModVector &p_mods,
                    ModVector &r_mods);

  // base class method. Adds linear inequalities
  void relaxInitFull(RelaxationPtr rel, bool *is_inf);

  // Does nothing.
  void relaxInitInc(RelaxationPtr rel, bool *is_inf);

  // Does nothing.
  void relaxNodeFull(NodePtr node, RelaxationPtr rel, bool *is_inf);

  // Does nothing.
  void relaxNodeInc(NodePtr node, RelaxationPtr rel, bool *is_inf);


  // base class method. Adds linearlization cuts when available.
  void separate(ConstSolutionPtr sol, NodePtr node, RelaxationPtr rel, 
                CutManager *cutman, SolutionPoolPtr s_pool, bool *sol_found,
                SeparationStatus *status);

  // base class method. 
  void writeStats(std::ostream &out) const;

private:
  /// Store statistics of presolving.
  struct SepaStats 
  {
    int iters;   ///> Number of times separation routine called. 
    int cuts;    ///> Number of cuts added.
    double time; ///> Total time used in separation
  };

  /// Store statistics of presolving.
  struct PresolveStats 
  {
    int iters;   ///> Number of iterations (main cycle).
    double time; ///> Total time used in initial presolve.
    double timeN;///> Total time used in presolveNode.
    int vBnd;    ///> Number of times variable-bounds were tightened.
    int nMods;   ///> Number of changes made in all nodes.
  };

  /// Absolute feasibility tolerance
  double aTol_;

  /// Logger.
  LoggerPtr logger_;
      
  /// For printing messages.
  static const std::string me_;

  /// Transformed problem (not the relaxation).
  ProblemPtr p_;

  /// Statistics about presolve
  PresolveStats pStats_;

  /// Relative feasibility tolerance
  double rTol_;

  /// Statistics about separation
  SepaStats sStats_;

  /// Keep track of time
  const Timer* timer_;

  /**
   * \brief Container for all bilinear functions. This should contain
   * variables and constraints of the problem, not the relaxation.
   */
  LinBilSet x0x1Funs_;

  /**
   * \brief Container for all square functions. This should contain variables and
   * constraints of the problem, not the relaxation.
   */
  LinSqrMap x2Funs_;

  /**
   * \brief Add a gradient-based linearization inequality.
   * \param[in] x       The variable x in (y = x^2)
   * \param[in] y       The variable y in (y = x^2)
   * \param[in] xl      x coordinate at which the gradient is evaluated
   * \param[in] yl      y coordinate at which the gradient is evaluated
   * \param[in] xval    x coordinate of point that is to be cut off
   * \param[in] yval    y coordinate of point that is to be cut off
   * \param[in] rel     Relaxation pointer to which the cut is added
   * \param[out] ifcuts True if the new inequality cuts off the point
   *                    (xval,yval)
   */
  void addCut_(VariablePtr x, VariablePtr y, double xl, double yl, double xval,
               double yval, RelaxationPtr rel, bool &ifcuts);

  /**
   * \brief Find the point at which a gradient-based linearization inequality
   * can be added.
   * \param[in] xval x coordinate of point that we want to cut off
   * \param[in] yval y coordinate of point that we want to cut off
   * \param[in] xl x coordinate of point at which gradient can be evaluated
   * \param[in] yl y coordinate of point at which gradient can be evaluated
   */
  void findLinPt_(double xval, double yval, double &xl, double &yl);

  /**
   * \brief Get one of the four linear functions and right hand sides for the
   * linear relaxation of a bilinear constraint y = x0x1.
   * \param[in] x0 x0 variable
   * \param[in] lb0 New lower bound of x0
   * \param[in] ub0 New upper bound of x0
   * \param[in] x1 x1 variable
   * \param[in] lb1 New lower bound of x1
   * \param[in] ub1 New upper bound of x1
   * \param[in] y The y variable
   * \param[in] type Could be 0,1,2,3. It tells the function which (out of the
   * four) linear constraints do we want. 
   * \param[out] rhs The rhs (upperbound) of the linear constraint
   * \return A linear function such that lf <= rhs is a relaxation of the
   * bilinear constraint.
   */
  LinearFunctionPtr getNewBilLf_(VariablePtr x0, double lb0, double ub0,
                                 VariablePtr x1, double lb1, double ub1,
                                 VariablePtr y, int type, double &rhs);

  /**
   * \brief Get linear function and right hand side for the linear
   * overestimator constraint for the square type y=x^2.
   * \param[in] x The variable x
   * \param[in] y The variable y
   * \param[in] lb The new lower bound of x
   * \param[in] ub The new upper bound of x
   * \param[out] r The rhs (upper bound) of the linear constraint.
   * \return A linear function such that lf <= rhs is a relaxation of the
   * square constraint.
   */
  LinearFunctionPtr getNewSqLf_(VariablePtr x, VariablePtr y, 
                                double lb, double ub, double & r);


  /// Return true if xval is one of the bounds of variable x
  bool isAtBnds_(ConstVariablePtr x, double xval);

  /**
   * \brief Strengthen bounds of variables in a bilinear constraint y=x0x1
   * \param[in] lx0x1 The bilinear term
   * \param[out] changed True if any bounds of y, x0 or x1 have changed. False
   * otherwise.
   * Return true if the new bounds are inconsistent (i.e., lb > ub for any of
   * the three variables)
   */
  bool propBilBnds_(LinBil* lx0x1, bool *changed);

  /**
   * \brief Strengthen bounds of variables in a bilinear constraint y=x0x1,
   * and save the modifications if any.
   * \param[in] lx0x1 The bilinear term
   * \param[in] rel The relaxation that is currently being solved.
   * \param[in] mod_rel If true, then change the relaxation also. The original
   * (or transformed) problem p_ is always modified.
   * \param[out] changed True if any bounds of y, x0 or x1 have changed. False
   * otherwise.
   * \param[in] p_mods A vector to save the modifications to the problem 
   * \param[in] r_mods A vector to save the modifications to the relaxation 
   * Return true if the new bounds are inconsistent (i.e., lb > ub for any of
   * the three variables)
   */
  bool propBilBnds_(LinBil* lx0x1, RelaxationPtr rel, bool mod_rel,
                    bool *changed, ModVector &p_mods, ModVector &r_mods);

  /**
   * \brief Strengthen bounds of variables in a square constraint y=x0^2.
   * \param[in] lx2 The sqaure term
   * \param[out] changed True if any bounds of y or x0 have changed. False
   * otherwise.
   * Return true if the new bounds are inconsistent (i.e., lb > ub for any of
   * the two variables)
   */
  bool propSqrBnds_(LinSqrMapIter lx2, bool *changed);

  /**
   * \brief Strengthen bounds of variables in a square constraint y=x0^2,
   * and save the modifications if any.
   * \param[in] lx2 The square term
   * \param[in] rel The relaxation that is currently being solved.
   * \param[in] mod_rel If true, then change the relaxation also. The original
   * (or transformed) problem p_ is always modified.
   * \param[out] changed True if any bounds of y or x0 have changed. False
   * otherwise.
   * \param[in] p_mods A vector to save the modifications to the problem 
   * \param[in] r_mods A vector to save the modifications to the relaxation 
   * Return true if the new bounds are inconsistent (i.e., lb > ub for any of
   * the two variables)
   */
  bool propSqrBnds_(LinSqrMapIter lx2, RelaxationPtr rel, bool mod_rel,
                    bool *changed, ModVector &p_mods, ModVector &r_mods);

  /**
   * \brief Relax all the square constraints: y=x0^2 and bilinear constraints
   * y=x0x1 in the problem.
   * \param[in] rel The relaxation to which the new linear inequalities are
   * added.
   * \param[out] is_inf True if the relaxation is detected to be infeasible
   * (because of inconsistent bounds or other reasons).
   */
  void relax_(RelaxationPtr rel, bool *is_inf);

  /// Reset all statistics to zero.
  void resetStats_();

  /**
   * \brief Modify bounds of a variable in the problem to the new bounds lb
   * and ub if the new bounds are tighter.
   * \param[in] v The variable
   * \param[in] lb The new lower bound
   * \param[in] ub The new upper bound
   * \param[out] changed True if the new bounds are tighter than the existing
   * bounds. False otherwise.
   * \return the number of bound changes (1 or 2) if any bounds are changed.
   * Returns -1 if the new bounds make the problem infeasible. Returns 0
   * otherwise.
   */
  int updatePBounds_(VariablePtr v, double lb, double ub, bool *changed);

  /**
   * \brief Modify bounds of a variable in the problem to the new bounds lb
   * and ub if the new bounds are tighter. Additionally, save the
   * modifications in the appropriate vectors so that they can be reverted at
   * a later time.
   * \param[in] v The variable
   * \param[in] lb The new lower bound
   * \param[in] ub The new upper bound
   * \param[in] rel The relaxation that is currently being solved
   * \param[in] mod_rel If true, then the relaxation is also modified. The
   * original (or transformed) problem is always modified.
   * \param[out] changed True if the new bounds are tighter than the existing
   * bounds. False otherwise.
   * \param[in] p_mods A vector to save the modifications to the problem 
   * \param[in] r_mods A vector to save the modifications to the relaxation 
   * \return the number of bound changes (1 or 2) if any bounds are changed.
   * Returns -1 if the new bounds make the problem infeasible. Returns 0
   * otherwise.
   */
  int updatePBounds_(VariablePtr v, double lb, double ub, RelaxationPtr rel,
                     bool mod_rel, bool *changed, ModVector &p_mods,
                     ModVector &r_mods);

  /**
   * \brief Update linear relaxation of the bilinear constraints after some
   * bounds have changed. The function checks whether each of the four
   * constraints are binding at the required extreme points of the box. If
   * not, the constraints are updated.
   * \param[in] lx0x1 The bilinear term
   * \param[in] rel The relaxation which contains the four linear constraints
   * \param[in] r_mods A vector into which the modifications are appended so
   * that they can be reverted at a later time.
   */
  void upBilCon_(LinBil* lx0x1, RelaxationPtr rel, ModVector &r_mods);


  /**
   * \brief Update linear relaxation of the square constraints (y=x^2) after
   * some bounds have changed. The function checks whether the upper bounding
   * constraint is binding at the required extreme points of the box. If not,
   * the constraint is updated.
   * \param[in] con The linear constraint that can be updated
   * \param[in] x The variable x. 
   * \param[in] y The variable y. 
   * \param[in] rel The relaxation which contains the linear constraint
   * \param[in] r_mods A vector into which the modifications are appended so
   * that they can be reverted at a later time.
   */
  void upSqCon_(ConstraintPtr con, VariablePtr x, VariablePtr y, 
                RelaxationPtr rel, ModVector &r_mods);

  /**
   * \brief Tighten the bounds on variables of the original (or transformed)
   * problem on the basis of the quadratic constraints y=x0^2 or y=x0x1.
   * Usually called in the initial presolve.
   * \param[out] changed True if any bounds are strengthened
   * \return true if some bounds become inconsistent, making the problem
   * infeasible. Flase otherwise.
   */
  bool varBndsFromCons_(bool *changed);
};

/// Shared pointer to QuadHandler.
typedef boost::shared_ptr<QuadHandler> QuadHandlerPtr;
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
