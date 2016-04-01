//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file TransPoly.h
 * \brief Declare class for reformulating a polynomial problem.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURTRANSPOLY_H
#define MINOTAURTRANSPOLY_H

#include "Transformer.h"

namespace Minotaur {

class MultilinearTermsHandler;
class MonomialFunction;
class YEqMonomial;
class YEqCG;
typedef boost::shared_ptr<MonomialFunction> MonomialFunPtr;
typedef boost::shared_ptr<MultilinearTermsHandler> MultilinearTermsHandlerPtr;


/**
 * \brief Base class for reformulating a general nonconvex mixed-integer
 * problem into one that has univariate constraints and multilinear constraints. 
 */
class TransPoly : public Transformer {
public:

  /// Default Constructor.
  TransPoly();

  /// Constructor.
  TransPoly(EnvPtr env, ConstProblemPtr p);

  /// Destroy.
  ~TransPoly();

  // base class method.
  std::string getName() const;

  // base class method.
  SolutionPtr getSolOrig(ConstSolutionPtr sol, int &err);

  // base class method.
  SolutionPtr getSolTrans(ConstSolutionPtr sol, int &err);

  // base class method.
  void reformulate(ProblemPtr &newp, HandlerVector &handlers, int &status);

private:

  /// Handler that takes care of constraints of the form y=u.v.w.x
  MultilinearTermsHandlerPtr mHandler_;

  /// Save monomials so that auxiliary variables can be reused.
  YEqMonomial *yMonoms_;

  /**
   * \brief Assign a constraint to multilinear-terms handler.
   * \param [in] con Constraint that is to be assigned to the handler. It must be
   * of the form y = x1.x2.x3...
   * \param [in] v The auxiliary variable y.
   * \param [in] mf The multilinear term saved as monomial. Care must be taken
   * to make sure all variables have power one in the monomial. See monomToMl_
   * function.
   */
  void assignMH_(ConstraintPtr con, VariablePtr v, MonomialFunPtr mf);

  /**
   * \brief Convert a monomial to another monomial with all powers equal to
   * one. Variables with higher powers get replaced by new variables and
   * constraints.
   * \param [in] mf The input monomial of the form
   * \f$x_1^{a_1}x_2^{a_2}\ldots\f$.
   * \return A monomial function that is of the form \f$y_1y_2\ldots\f$.
   */
  MonomialFunPtr monomToMl_(MonomialFunPtr mf);

  /**
   * \brief Assign an auxiliary variable for a monomial.
   *
   * If the monomial was
   * already handled before, return the associated variable, otherwise new
   * constraint and variable are added. If the monomial has powers higher than
   * one on some variables, additional transformation is carried out.
   * \param [in] mf The input monomial of the form
   * \f$x_1^{a_1}x_2^{a_2}\ldots\f$.
   * \return A variable \f$y\f$ corresponding to
   * \f$y = x_1^{a_1}x_2^{a_2}\ldots\f$
   */
  VariablePtr newPolyVar_(MonomialFunPtr mf);

  /**
   * \brief Assign an auxiliary variable for a monomial, linear function or a
   * product of variable and constant.
   *
   * This function may call newPolyVar_ above. This function judiciously
   * creates auxiliary variables for the input monomial or linear functions.
   * \param [in] mf The input monomial of the form
   * \f$x_1^{a_1}x_2^{a_2}\ldots\f$. It may be NULL.
   * \param [in] lf The input linear function. It may be NULL.
   * \param [in] v The input variable. We want to have a constraint
   * \f$y = kv\f$ in this case. mf, lf, v should be all NULL except one.
   * \param [in] d double value that must be added to the mf or lf or the
   * variable.
   * \param [in] k The multiplier of v.
   * \return The new auxiliary variable.
   */
  VariablePtr newPolyVar_(const CNode *cnode, MonomialFunPtr mf,
                          LinearFunctionPtr lf, VariablePtr v,
                          double d, double k);

  /**
   * \brief Traverse the computational graph recursively and transform the
   * nonlinear constraint.
   *
   * \param [in] node The node whose sub-tree must be reformulated.
   * \param [out] mf The monomial function representing this node.
   * \param [out] v The variable representing this node.
   * mf, lf, v should be all NULL except one.
   * \param [in] d double value that is added to the mf, lf, or v at this node.
   * \param [in] k The multiplier of v.
   */
  void recursPolyRef_(const CNode *node, 
                      MonomialFunPtr &mf, LinearFunctionPtr &lf,
                      VariablePtr &v, double &d, double &k);

  /**
   * \brief Reformulate the node with OpMinus operation in the computational
   * graph.
   *
   * Depending on the type of functions in the two children, this function
   * adds new variables, or combines linear functions for the reformulation.
   * \param [in] mfl Monomial function from the left child. It may be NULL.
   * \param [in] mfr Monomial function from the right child. It may be NULL.
   * \param [in] lfl Linear function from the left child. It may be NULL
   * \param [in] lfr Linear function from the right child. It may be NULL
   * \param [in] vl Variable from the left child.
   * \param [in] vr Variable from the right child.
   * \param [in] dl double value from the left child.
   * \param [in] dr double value from the right child.
   * \param [in] kl The multiplier of vl.
   * \param [in] kr The multiplier of vr.
   * \param [out] mf The new monomial function. It may be NULL.
   * \param [out] lf The new linear function. It may be NULL.
   * \param [out] v The new variable. It may be NULL.
   * \param [out] d The new double value for addition.
   * \param [out] k The multiplier of v.
   */
  void refMinus_(MonomialFunPtr mfl, MonomialFunPtr mfr,
                 LinearFunctionPtr lfl, LinearFunctionPtr lfr,
                 VariablePtr vl, VariablePtr vr,
                 double dl, double dr,
                 double kl, double kr,
                 MonomialFunPtr &mf, LinearFunctionPtr &lf,
                 VariablePtr &v, double &d, double &k);

  /**
   * \brief Reformulate the node with OpMult operation in the computational
   * graph.
   *
   * Depending on the type of functions in the two children, this function
   * adds new variables, or combines linear functions for the reformulation.
   * \param [in] mfl Monomial function from the left child. It may be NULL.
   * \param [in] mfr Monomial function from the right child. It may be NULL.
   * \param [in] lfl Linear function from the left child. It may be NULL
   * \param [in] lfr Linear function from the right child. It may be NULL
   * \param [in] vl Variable from the left child.
   * \param [in] vr Variable from the right child.
   * \param [in] dl double value from the left child.
   * \param [in] dr double value from the right child.
   * \param [in] kl The multiplier of vl.
   * \param [in] kr The multiplier of vr.
   * \param [out] mf The new monomial function. It may be NULL.
   * \param [out] lf The new linear function. It may be NULL.
   * \param [out] v The new variable. It may be NULL.
   * \param [out] d The new double value for addition.
   * \param [out] k The multiplier of v.
   */
  void refMult_(MonomialFunPtr mfl, MonomialFunPtr mfr,
                LinearFunctionPtr lfl, LinearFunctionPtr lfr,
                VariablePtr vl, VariablePtr vr,
                double dl, double dr,
                double kl, double kr,
                MonomialFunPtr &mf, LinearFunctionPtr &lf,
                VariablePtr &v, double &d, double &k);

  
  /*
   * \brief Reformulate the constraints of the original problem into new
   * problem.
   */
  void refNonlinCons_();
 
  /*
   * \brief Reformulate the objective function of the original problem into new
   * problem.
   */
  void refNonlinObj_();
    
  /**
   * \brief Reformulate the node with OpPlus operation in the computational
   * graph.
   *
   * Depending on the type of functions in the two children, this function
   * adds new variables, or combines linear functions for the reformulation.
   * \param [in] mfl Monomial function from the left child. It may be NULL.
   * \param [in] mfr Monomial function from the right child. It may be NULL.
   * \param [in] lfl Linear function from the left child. It may be NULL
   * \param [in] lfr Linear function from the right child. It may be NULL
   * \param [in] vl Variable from the left child.
   * \param [in] vr Variable from the right child.
   * \param [in] dl double value from the left child.
   * \param [in] dr double value from the right child.
   * \param [in] kl The multiplier of vl.
   * \param [in] kr The multiplier of vr.
   * \param [out] mf The new monomial function. It may be NULL.
   * \param [out] lf The new linear function. It may be NULL.
   * \param [out] v The new variable. It may be NULL.
   * \param [out] d The new double value for addition.
   * \param [out] k The multiplier of v.
   */
  void refPlus_(MonomialFunPtr mfl, MonomialFunPtr mfr,
                LinearFunctionPtr lfl, LinearFunctionPtr lfr,
                VariablePtr vl, VariablePtr vr,
                double dl, double dr,
                double kl, double kr,
                MonomialFunPtr &mf, LinearFunctionPtr &lf,
                VariablePtr &v, double &d, double &k);

  /**
   * \brief Reformulate univariate nonlinear functions like square-root,
   * constant power, log, sin etc.
   *
   * This function is called anytime a nonlinear operation is seen in the
   * traversal of computational graph. If the expression is constant, it
   * complains. Otherwise, newPolyVar_ function is called to reformulate.
   * \param [in] node The node whose sub-tree must be reformulated.
   * \param [in] mfl The monomial function representing this node.
   * \param [in] lfl The linear function representing this node.
   * \param [in] vl The variable representing this node.
   * \param [in] dl The addition-constant associated with this node.
   * \param [in] kl The multiplier of vl associated with this node.
   * \param [out] v Auxiliary variable for this node. May be NULL.
   * \param [out] d The addition-constant, if the sub-function is constant.
   * \param [out] k Multiplier of auxiliary variable v.
   */
  void refUnivarOpPoly_(const CNode *node, MonomialFunPtr mfl,
                        LinearFunctionPtr lfl, VariablePtr vl, double dl,
                        double kl, VariablePtr &v, double &d,
                        double &k);
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
