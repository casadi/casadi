//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file Constraint.h
 * \brief Get information about a constraint in a given Problem.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURCONSTRAINT_H
#define MINOTAURCONSTRAINT_H

#include <string>

#include "Types.h"

namespace Minotaur {
  class   Function;
  class   LinearFunction;
  class   QuadraticFunction;
  class   NonlinearFunction;
  typedef boost::shared_ptr<Function> FunctionPtr;
  typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;
  typedef boost::shared_ptr<NonlinearFunction> NonlinearFunctionPtr;
  typedef boost::shared_ptr<QuadraticFunction> QuadraticFunctionPtr;
  typedef std::set<std::pair<VariablePtr, FunctionType> >::const_iterator 
    VariableFunIterator;

  /**
   * \brief The Constraint class is used to manage a constraint.
   *
   * The general form of a constraint is:
   * \f[
   * l \leq f(x) \leq u,
   * \f]
   * where \f$f(x)\f$ is a function of the variables and \f$l, u\f$ are
   * constraints placed on \f$f\f$. \f$l\f$ can be \f$-\infty\f$ and/or
   * \f$u\f$ can be \f$\infty\f$.
   * 
   * A constraint should only be constructed through a Problem object. Doing
   * this ensures that constraints have ID's unique for that problem. A
   * constraint can thus be part of only one problem. A function \f$f\f$ that
   * is part of a constraint should also not be shared at several places in
   * the same problem or amongst different problems.  Similarly, constraints
   * should be deleted only through a problem.  Changes in a constraint (e.g.
   * changing name, bounds etc) may be done directly.
   * 
   * The state of a constraint describes whether the constraint is free or
   * deleted etc. It is Fixed if each variable in the constraints has been 
   * fixed.
   * 
   * A Constraint also stores a map of variables that occur in that
   * constraint. The map tells whether the variable is in a non-linear, 
   * or quadratic etc part of the constraint. A variable that occurs 
   * once in linear and thrice in quadratic parts will have four entries in
   * the map, once for each occurrence.
   * 
   * The name of a constraint can be set by the user. There are no checks for
   * duplicate names. The ID of the constraint however is unique in the
   * problem.
   */
  class Constraint {
    public:

      /// Only Problem class can modify a Constraint. All modification methods
      /// are private.
      friend class Problem;


      /// Default constructor.
      Constraint();

      /**
       * \brief This is the only real constructor. All others should call this
       * constructor.
       *
       * Create a constraint with 
       * \param [in] id the unique id of the constraint (determined by
       * Problem object), 
       * \param [in] index the index of the constraint in the problem,
       * \param [in] f the function in the body of the constraint,
       * \param [in] lb the lower bound, can be -INFINITY,
       * \param [in] ub the upper bound, can be  INFINITY,
       * \param [in] name The string name of the constraint.
       */
      Constraint(UInt id, UInt index, FunctionPtr f, double lb, double ub, 
                 std::string name);

      /// Destroy
      virtual ~Constraint();

      /// Get the value or activity at a given point.
      double getActivity(const double *x, int *error) const;

      /// Return a pointer to the function.
      const FunctionPtr getFunction() const { return f_; }

      /// Get the function type.
      FunctionType getFunctionType() const;

      /// Return the unique id of the constraint.
      UInt getId() const { return id_; }

      /// Return the index of the constraint.
      UInt getIndex() const { return index_; }

      /// Get the 'l' value. or the lower bound constraint on 'f'.
      double getLb() const { return lb_; }

      /// Get the linear part of the constraint function 'f'.
      const LinearFunctionPtr getLinearFunction() const;

      // Get the name of the constraint. 
      const std::string getName() const;

      /// Get the nonlinear part of the constraint function 'f'.
      const NonlinearFunctionPtr getNonlinearFunction() const;

      /// Get the quadratic part of the constraint function 'f'.
      const QuadraticFunctionPtr getQuadraticFunction() const;

      /// Get the current state of the constraint: freed, fixed etc.
      ConsState getState() const { return state_; }

      /// Get the 'u' value. or the upper bound constraint on 'f'.
      double getUb() const { return ub_; }

      /// display the constraint
      void write(std::ostream &out) const;

    protected:
      /// Add a constant to the constraint function. lb - c <= f <= ub - c.
      void add_(double cb);

      /// Change the linear part of constraint.
      void changeLf_(LinearFunctionPtr lf);

      /// Change the nonlinear part of constraint.
      void changeNlf_(NonlinearFunctionPtr nlf);

      /// Delete variables fixed at value val.
      void delFixedVar_(VariablePtr v, double val);

      /**
       * \brief Negate the constraint.
       *
       * Convert the constraint to 
       * \f[
       * -u \leq -f(x) \leq -l
       *  \f]
       */
      void reverseSense_();

      /// Set the id. 
      void setId_(UInt n) { id_ = n; }

      /// Set the index. 
      void setIndex_(UInt n) { index_ = n; }

      /**
       * \brief Set a new lower bound.
       *
       * The new lower bound can be -INFINITY. It can change the state of
       * the constraint.
       */
      void setLb_(double newlb) { lb_ = newlb; }

      /// Set name of the constraint
      void setName_(std::string name);

      /// Set state of the constraint.
      void setState_(ConsState state) { state_ = state; return; }

      /**
       * \brief Set a new upper bound.
       *
       * The new bound can be INFINITY. It can affect the state of the
       * constraint.
       */
      void setUb_(double newub) { ub_ = newub; }

      /// \brief Substitute a variable \f$x_1\f$ by \f$rx_2\f$.
      void subst_(VariablePtr out, VariablePtr in, double rat, bool *instay);

    private:
      /// The function 'f' in l <= f(x) <= u.
      FunctionPtr f_;

      /// id that is unique for this constraint in the problem.
      UInt id_;         

      /// id that is unique for this constraint in the problem.
      UInt index_;         

      /// 'l' [-infinity, infinity).
      double lb_;

      /// name of the constraint. could be NULL.
      std::string name_;

      /// free or fixed etc.
      ConsState state_;

      /// 'u' (-infinity, infinity].
      double ub_;
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
