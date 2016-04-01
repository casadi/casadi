// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file AMPLNonlinearFunction.h
 * \brief Declare the AMPLNonlinearFunction class for setting up evaluation
 * and derivatives routines of nonlinear functions through AMPL.
 * functions from AMPL.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURAMPLNONLINEARFUNCTION_H
#define MINOTAURAMPLNONLINEARFUNCTION_H

#include "asl.h"

#include "NonlinearFunction.h"

namespace MINOTAUR_AMPL {

/**
 * \brief Declare the AMPLNonlinearFunction class for setting up evaluation
 * and derivatives of nonlinear Functions.
 *
 * This class does not contain the computational graph of the function. Use it
 * only for evaluating function values and derivatives.
 */
class AMPLNonlinearFunction : public Minotaur::NonlinearFunction {
public:

  /// Default constructor.
  AMPLNonlinearFunction ();

  /**
   * \brief Create a nonlinear function for a given constraint.
   * \param [in] i The constraint number of which we need the nonlinear
   * function. This value is ignored if is_obj is True.
   * \param [in] nvars Total number of variables in the problem.
   * \param [in] my_asl Pointer to ASL for calling its routines.
   * \param [in] is_obj If True, then get function from objective. Otherwise,
   * get function of i-th constriant.
   */
  AMPLNonlinearFunction (Minotaur::UInt i, Minotaur::UInt nvars, ASL* my_asl,
                     bool is_in_obj);

  // Not available.
  Minotaur::NonlinearFunctionPtr 
    cloneWithVars(Minotaur::VariableConstIterator vbeg, 
                  int *err) const;

  // Evaluate at a point. Base class function.
  double eval(const double *x, int *error);

  // Evaluate gradient at a point. Base class function.
  void evalGradient(const double *x, double *grad_f,
                    int *error);

  void evalHessian(const double mult, const double *x, 
                   const Minotaur::LTHessStor *stor, double *values, 
                   int *error);

  // Not available.
  void  fillHessStor(Minotaur::LTHessStor *);

  // Not available.
  void  finalHessStor(const Minotaur::LTHessStor *);

  // Not available.
  void fillJac(const double *, double *, int *);

  // Get variables used in this function. Base class method.
  void getVars(Minotaur::VariableSet *);

  // Multiply by a constant. Base class method.
  void multiply(double c);

  // Not available.
  void prepJac(Minotaur::VarSetConstIter, Minotaur::VarSetConstIter);

  /**
   * \brief Tell what variables are in this function.
   *
   * These variabvles are then
   * stored in this class for future use by Minotaur routines.
   *
   * \param [in] vb Iterator pointing to first element of set.
   * \param [in] ve Iterator pointing to end of set.
   */
  void setVars(Minotaur::VarSetConstIterator vb, 
               Minotaur::VarSetConstIterator ve);

  // Display the function. Base class method.
  void write(std::ostream &out) const;

private:
  /// Number of variables in the problem;
  Minotaur::UInt nVars_; 

  /// Index of the corresponding constraint/obj in AMPL.
  Minotaur::UInt amplIndex_; 

  // pointer to ampl's asl.
  ASL *myAsl_;

  /**
   * \brief True, if the function is negated, e.g. when maximizing instead of
   * minimizing.
   */
  bool neg_;              

  /// True if it appears in objective.
  bool isInObj_;
}; // AMPLNonlinearFunction

typedef boost::shared_ptr<AMPLNonlinearFunction> AMPLNlfPtr;
}  // namespace

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
