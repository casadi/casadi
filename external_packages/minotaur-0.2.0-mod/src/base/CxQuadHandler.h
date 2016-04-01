//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file CxQuadHandler.h
 * \brief Define the CxQuadHandler class for handling convex quadratic
 * objective functions and constraints.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURCXQUADRATICHANDLER_H
#define MINOTAURCXQUADRATICHANDLER_H

#include "Handler.h"

namespace Minotaur {

class Engine;
class Function;
class LinearFunction;
class Objective;
class Problem;
class QuadraticFunction;
typedef boost::shared_ptr<Engine> EnginePtr;
typedef boost::shared_ptr<Function> FunctionPtr;
typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;
typedef boost::shared_ptr<Objective> ObjectivePtr;
typedef boost::shared_ptr<const Problem> ConstProblemPtr;
typedef boost::shared_ptr<QuadraticFunction> QuadraticFunctionPtr;

/// Save information about constraints of the form \f$ y \leq x^2 \f$.
struct Secant {
  VariablePtr   auxVar;   /// The variable y.
  VariablePtr   sqVar;    /// The variable x.
  ConstraintPtr cons;     /// The linear constraint that is used as approx.
};

/// Pointer to Secant
typedef Secant * SecantPtr;

/// Vector-iterator for Secant
typedef std::vector< SecantPtr >::iterator SecantIterator;

/// Map of 'x' and the secant that is used for \f$ y \leq x^2 \f$.
typedef std::map<VariablePtr, SecantPtr, CompareVariablePtr> VarSecantMap;

/// Iterator for VarSecantMap
typedef VarSecantMap::iterator VarSecantMapIter;


/**
 * A McCormick object stores some information about McCormick inequalities 
 * for the constraints of the form \f$x_0x_1 \leq y\f$ or \f$x_0x_1 \geq y \f$
 * or \f$x_0x_1 = y\f$.
 */
class McCormick {
public:
  /// LT: x0x1 <= y; GT: x0x1 >= y; EQ: x0x1 = y
  typedef enum {
    LT,
    GT,
    EQ
  } Sense;

private:
  /// First variable.
  VariablePtr x0_;

  /// Second variable.
  VariablePtr x1_;

  /// Auxiliary variable.
  VariablePtr y_;

  /// Greater than, less than or equal to.
  Sense s_;

  /// Constraint 0.
  ConstraintPtr c0_;
        
  /// Constraint 1.
  ConstraintPtr c1_;

  /// Constraint 2.
  ConstraintPtr c2_;

  /// Constraint 3.
  ConstraintPtr c3_;

public:
  /// Default constructor. 
  McCormick(VariablePtr x0, VariablePtr x1, Sense sense);

  /// Destroy.
  ~McCormick();

  /**
   * Set the auxiliary variable. If a bilinear term repeats itself, all
   * such terms must share the same auxiliary variable.
   */
  void setAux(VariablePtr y) {y_ = y;};

  /// Set the zeroth constraint: \f$ y \geq l_0x_1 + l_1x_0 - l_1l_0 \f$.
  void setC0(ConstraintPtr c) {c0_ = c;};

  /// Set the first constraint: \f$ y \geq u_0x_1 + u_1x_0 - u_1u_0 \f$.
  void setC1(ConstraintPtr c) {c1_ = c;};

  /// Set the third constrait: \f$ y \leq u_1x_0 + l_0x_1 - l_0u_1 \f$.
  void setC2(ConstraintPtr c) {c2_ = c;};

  /// Set the third constrait: \f$ y \leq l_1x_0 + u_0x_1 - l_1u_0 \f$.
  void setC3(ConstraintPtr c) {c3_ = c;};

  /// Get one of the four constraints.
  ConstraintPtr getC0() {return c0_;};
  ConstraintPtr getC1() {return c1_;};
  ConstraintPtr getC2() {return c2_;};
  ConstraintPtr getC3() {return c3_;};

  /// Get the auxiliary variable.
  VariablePtr getAux() {return y_;};

  /// Get \f$x_0\f$
  VariablePtr getX0() {return x0_;};

  /// Get \f$x_1\f$
  VariablePtr getX1() {return x1_;};

  /// Get the variable other than x, in the product.
  VariablePtr getOtherX(ConstVariablePtr x) const;

  /// Get the sense of the bilinear constraint: LT, EQ, GT.
  Sense getSense() {return s_;};

  /// Set the sense of the bilinear constraint: LT, EQ, GT.
  void  setSense(Sense sense) {s_ = sense;};

  /// Check if a bilinear constraint is violated at the current point x.
  bool isViolated(const double *x, const double &tol) const;

  /**
   * Check if a bilinear constraint is violated for the given values of
   * \f$x_0, x_1, y\f$.
   */
  bool isViolated(const double &x0val, const double &x1val, 
                  const double &y0val, const double &tol) const;
};
/// shared pointer to McCormick object.
typedef boost::shared_ptr<McCormick> McCormickPtr;

/**
 * Compare two McCormick objects. Since we keep them in a set, we need to
 * sort them. We use lexicographic ordering (i.e. based on ids of 
 * \f$(x_0, x_1)\f$).
 */
struct CompareMcCormick {
  bool operator()(McCormickPtr b0, McCormickPtr b1) const;
};

/// A set of McCormick objects.
typedef std::set<McCormickPtr, CompareMcCormick> McCormickSet;

/// Iterator of McCormick objects over a set.
typedef McCormickSet::iterator McCormickSetIter;


/**
 * An CxQuadHandler handles the convex parts of quadratic functions of a
 * problem. For now, we will just handle squares of singleton variables e.g.
 * \f$\sum_ix_i^2 \leq u_0\f$. \f$u_0\f$ could be an auxiliary varible or a
 * constant. Later we will introduce sums of squares of
 * general linear functions as well.
 */
class CxQuadHandler : public Handler {
protected:
  /**
   *  For each constraint of the type \f$y \leq x^2\f$, we add a new
   *  secant approximation. This map stores the 'y' variables and the
   *  associated linear secant-constraint.
   */
  VarSecantMap cvCons_;

  /**
   * Keep a set of McCormick inequalities that were added and for which
   * we will update the co-efficients etc.
   */
  McCormickSet mcCons_;

  /**
   * Variables that occur in bilinear terms and also concave square terms. 
   * These do not include auxiliary variables that are added in relaxation.
   * These are all the candidates that we consider for branching.
   */
  VarSet brVars_;

  /// Tolerance.
  double eTol_;

  /// Logger.
  LoggerPtr logger_;
      
  /// For printing.
  static const std::string me_;

  /// Original problem.
  ProblemPtr problem_;

  /**
   * Add quadratic/linear relaxations of the quadratic range constraint 
   * 'cons'.
   */
  void relaxTwoSided_(QuadraticFunctionPtr qf, ConstraintPtr cons,
                      RelaxationPtr rel);

  /**
   * Add quadratic/linear relaxations of the quadratic constraint 'cons'
   * that only has an upperbound.
   */
  void relaxOneSided_(QuadraticFunctionPtr qf, ConstraintPtr cons, 
                      RelaxationPtr rel);

  /// Relax the objective function, to make it convex 
  void relaxObj_(ObjectivePtr obj, RelaxationPtr rel);

  /**
   * Get secant approximation for the inequality:
   * \f$ y - x^2 \leq 0\f$, and add it to the relaxation.
   */
  void addSecant_(VariablePtr x, VariablePtr y, RelaxationPtr rel);

  /// Get linear function and right hand side (r) for a secant constraint.
  LinearFunctionPtr getNewSecantLf_(VariablePtr x, VariablePtr y, 
                                    double & lb, double & ub, double & r);

  /// Add all four McCormick inequalities for \f$ y = x_0x_1\f$.
  VariablePtr addMcCormick_(VariablePtr x0, VariablePtr x1, 
                            RelaxationPtr rel);

  /// Add two McCormick inequalities for \f$ y \geq x_0x_1\f$.
  VariablePtr addMcCormickLower_(VariablePtr x0, VariablePtr x1, 
                                 RelaxationPtr rel);

  /// Add two McCormick inequalities for \f$ y \leq x_0x_1\f$.
  VariablePtr addMcCormickUpper_(VariablePtr x0, VariablePtr x1, 
                                 RelaxationPtr rel);

  /// Generate the appropriate McCormick inequality using the bounds.
  LinearFunctionPtr getMcLf_(VariablePtr x0, double lb0, double ub0,
                             VariablePtr x1, double lb1, double ub1, VariablePtr y, 
                             double &rhs, UInt i);

  void binToLin_();
  void binToLinFun_(FunctionPtr f, LinearFunctionPtr lf2);

  /**
   * For now, this handler will introduce a new variable for each sum of
   * squares of variables. We will also add constraints, that have these
   * functions, from the original problem if such constraints have not
   * already been added. If these constraints are already in the problem,
   * we will add the corresponding new variable in the linear function of
   * the constraint. Bounds: \f$ u_O \in [0,\sum_i\max(lb_i^2, ub_i^2)\f$].
   */
  void relax_(RelaxationPtr rel, bool *is_inf);

  void removeFixed_();
  void removeFixedFun_(FunctionPtr f, LinearFunctionPtr lf2, double *c);
public:
  /// Default constructor.
  CxQuadHandler(EnvPtr env, ProblemPtr problem);

  /// Destroy
  ~CxQuadHandler();

  // Does nothing.
  void relaxInitFull(RelaxationPtr rel, bool *is_inf);

  // Does nothing.
  void relaxInitInc(RelaxationPtr rel, bool *is_inf);

  // Does nothing.
  void relaxNodeFull(NodePtr node, RelaxationPtr rel, bool *is_inf);

  // Does nothing.
  void relaxNodeInc(NodePtr node, RelaxationPtr rel, bool *is_inf);


  /**
   * Suppose we added a variable \f$u_0 \geq \sum_i x_i^2\f$. In this
   * function, check if at a given point, \f$(u_0^*, x^*)\f$, the above
   * constraint holds or not. Return true if the constraint holds, false
   * otherwise.  Checks the conditions for all such constraints but stops
   * at the first infeasible one.
   */
  bool isFeasible(ConstSolutionPtr sol, RelaxationPtr relaxation, 
                  bool &is_inf, double &inf_meas);

  /**
   * Not implemented yet.
   */
  void separate(ConstSolutionPtr sol, NodePtr node, RelaxationPtr rel, 
                CutManager *cutman, SolutionPoolPtr s_pool, bool *sol_found,
                SeparationStatus *status);


  /// Return \f$u_0\f$ it is constrained to be an integer.
  void getBranchingCandidates(RelaxationPtr rel, 
                              const DoubleVector &x, ModVector &mods,
                              BrVarCandSet &cands, BrCandVector &,
                              bool & is_inf);

  // Implement Handler::getBrMod().
  ModificationPtr getBrMod(BrCandPtr cand, DoubleVector &x, 
                           RelaxationPtr rel, BranchDirection dir);

  // Implement Handler::getBranches().
  Branches getBranches(BrCandPtr cand, DoubleVector & x,
                       RelaxationPtr rel, SolutionPoolPtr s_pool);

  // presolve.
  SolveStatus presolve(PreModQ *pre_mods, bool *changed);

  // Implement Handler::presolveNode().
  bool presolveNode(RelaxationPtr rel, NodePtr node,
                    SolutionPoolPtr s_pool, ModVector &p_mods,
                    ModVector &r_mods);

  // Write name
  std::string getName() const;

};

/// Shared pointer to CxQuadHandler.
typedef boost::shared_ptr<CxQuadHandler> CxQuadHandlerPtr;

/// Shared pointer to const CxQuadHandler.
typedef boost::shared_ptr<const CxQuadHandler> CxQuadConstHandlerPtr;
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
