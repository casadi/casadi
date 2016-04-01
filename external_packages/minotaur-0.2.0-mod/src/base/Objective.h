// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file Objective.h
 * \brief Declare the Objective class for storing and manipulating objective
 * function.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef OBJECTIVE_H
#define OBJECTIVE_H

#include <string>

#include "Types.h"

namespace Minotaur {

  class   Function;
  class   LinearFunction;
  class   QuadraticFunction;
  class   NonlinearFunction;
  typedef boost::shared_ptr<Function> FunctionPtr;
  typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;
  typedef boost::shared_ptr<const LinearFunction> ConstLinearFunctionPtr;
  typedef boost::shared_ptr<NonlinearFunction> NonlinearFunctionPtr;
  typedef boost::shared_ptr<QuadraticFunction> QuadraticFunctionPtr;
  typedef std::set<std::pair<VariablePtr, FunctionType> >::const_iterator 
    VariableFunIterator;
  /**
   * The Objective class is meant to implement an objective function for the
   * problem (or formulation). There are only a few differences between
   * Objective and Constraint. The Objective has a sense (Maximize or
   * Minimize) associated with it. The objective may also have a constant
   * term. Like a constraint, an objective has a function and lower and upper
   * bounds. These bounds are optional. An objective must have the ability to
   * be evaluated at a given point, and also to return its gradient at a given
   * point.
   * 
   * Minotaur will always minimize an objective internally. So if the
   * objective is to maximize, it will try to minimize the negative.
   */
  class Objective
  {
    public:
      /// Only Problem class can modify an Objective. All modification methods
      /// are private.
      friend class Problem;

      /// Default constructor
      Objective();

      /**
       * Create an objective with a specific id, constant term and sense. Does
       * not need a function. A default name is created.
       */
      Objective(double cb = 0.0, ObjectiveType otyp = Minimize);

      /*
       * Create an objective with a specific function, id, constant term and
       * sense. A default name is created.
       */
      Objective(FunctionPtr fPtr, double cb = 0.0, 
                ObjectiveType otyp = Minimize);

      /**
       * Create an objective with a specific function, id, constant term,
       * sense and a name.
       */
      Objective(FunctionPtr fPtr, double cb, ObjectiveType otyp, 
                std::string name);

      /// Destroy
      virtual ~Objective() { }


      // Functions that query or get information.

      /// Return the constant term in this objective.
      double getConstant() const { return cb_; }

      /// Return the type (Maximize or Minimize) of this objective.
      ObjectiveType getObjectiveType() const;

      /// Return the type of the function in objective.
      FunctionType getFunctionType() const;

      /// Return the state of the objective. Fixed, removed, freed etc.
      ObjState getState() const { return state_; }

      /*
       * Return the name of this objective. If the user did not create a name
       * for the objective, it is generated automatically
       */
      virtual const std::string getName() const;

      /// Return the function in the objective
      const FunctionPtr getFunction() const {
        return f_;
      }

      /**
       * Return the linear function part of the objective function. Could be
       * NULL.
       */
      const LinearFunctionPtr getLinearFunction() const;

      /**
       * Return the quadratic function part of the objective function. Could be
       * NULL.
       */
      const QuadraticFunctionPtr getQuadraticFunction() const;

      /**
       * Return the nonlinear function part of the objective function. Could be
       * NULL.
       */
      const NonlinearFunctionPtr getNonlinearFunction() const;

      /// Print the objective function to the output stream. 
      void write(std::ostream &out) const; 

      /**
       * Evaluate the objective function (along with the constant term) at the
       * given point x.
       */
      double eval(const double *x, int *err) const;

      /**
       * Evaluate the gradient at the given point x and fill in the gradient
       * values in the array grad_f. The array grad_f is assumed to have the
       * required size (equal to number of variables in the problem). 
       */
      void evalGradient(const double *x, double *grad_f, int *error);

      /**
       * Evaluate the gradient at the given point and fill in the gradient
       * values in the array grad_f. The array grad_f is assumed to have the
       * required size (equal to number of variables in the problem). 
       * \todo XXX: replace this function by the one below.
       */
      void evalGradient(VariableGroup &point, double *grad_f);

      // Functions that modify or change the objective.
      
    private:
      /**
       * Constant term in the objective. When the objective is evaluated at a
       * certain point, this term is included.
       */
      double cb_;

      /// Maximize or Minimize.
      ObjectiveType otyp_;

      /// The state of the function in the objective. Freed, fixed etc.
      ObjState state_;

      /// The function in the objective. Could be NULL.
      FunctionPtr f_;

      /**
       * The name of this objective. The user can set or change this. The
       * default name is derived from the id.
       */
      std::string name_;

      void delFixedVar_(VariablePtr v, double val);

      /// Set the state of the objective.
      void setState_(ObjState state) { state_ = state; return; }

      /** 
       * Take the sum of the existing linear function in the objective and
       * lPtr and set this new function as the linear function of objective.
       * Quadratic and other nonlinear parts are unaffected. lPtr is also
       * unaffected.
       */
      void add_(ConstLinearFunctionPtr lPtr);

      /// Add a constant to the objective function.
      void add_(double cb);

      /// Remove the quadratic part of the function and return a pointer to it.
      QuadraticFunctionPtr removeQuadratic_();

      void subst_(VariablePtr out, VariablePtr in, double rat=1.0);

      /**
       * Change the sense of objective and also multiply by (-1). Instead of
       * max f(x), the objective is now min -f(x).
       */
      void negate_();
  };

  typedef boost::shared_ptr<Objective> ObjectivePtr;
  typedef boost::shared_ptr<const Objective> ConstObjPtr;
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
