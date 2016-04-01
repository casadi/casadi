//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file LinBil.h
 * \brief Declare a class for storing linear under and overestimators of
 * bilinear functions
 * \f$ y_1 = x_1x_2 \f$,
 * \author Ashutosh Mahajan, IIT Bombay
 */

#ifndef MINOTAURLINBIL_H
#define MINOTAURLINBIL_H

#include "Handler.h"

namespace Minotaur {

/**
 * A LinBil object stores some information about linear-relaxation inequalities 
 * for the bilinear constraints of the form \f$x_0x_1 = y\f$ 
 */
class LinBil {
private:
  /// Absolute feasibility tolerance
  double aTol_;

  /// Constraint 0.
  ConstraintPtr c0_;
        
  /// Constraint 1.
  ConstraintPtr c1_;

  /// Constraint 2.
  ConstraintPtr c2_;

  /// Constraint 3.
  ConstraintPtr c3_;

  /// Relative feasibility tolerance
  double rTol_;

  /// First variable.
  VariablePtr x0_;

  /// Second variable.
  VariablePtr x1_;

  /// Auxiliary variable.
  VariablePtr y_;


public:
  /// Default constructor. 
  LinBil(VariablePtr y, VariablePtr x0, VariablePtr x1);

  /// Destroy.
  ~LinBil();

  /// Get the first out of the four constraints.
  ConstraintPtr getC0() {return c0_;};

  /// Get the second out of the four constraints.
  ConstraintPtr getC1() {return c1_;};

  /// Get the third out of the four constraints.
  ConstraintPtr getC2() {return c2_;};

  /// Get the fourth out of the four constraints.
  ConstraintPtr getC3() {return c3_;};

  /// Get the auxiliary variable.
  VariablePtr getY() {return y_;};

  /// Get \f$x_0\f$
  VariablePtr getX0() {return x0_;};

  /// Get \f$x_1\f$
  VariablePtr getX1() {return x1_;};

  /// Get the variable other than x, in the product.
  VariablePtr getOtherX(ConstVariablePtr x) const;

  /// Check if a bilinear constraint is violated at the current point x.
  bool isViolated(const double *x, double &vio) const;

  /**
   * \brief Check if a bilinear constraint is violated for the given values of
   * \f$x_0, x_1, y\f$.
   */
  bool isViolated(const double x0val, const double x1val, 
                  const double y0val) const;

  void setCons(ConstraintPtr c0, ConstraintPtr c1, ConstraintPtr c2,
               ConstraintPtr c3);
};

/**
 * Compare two LinBil objects. Since we keep them in a set, we need to
 * sort them. We use lexicographic ordering (i.e. based on ids of 
 * \f$(x_0, x_1)\f$).
 */
struct CompareLinBil {
  bool operator() (LinBil* b0, LinBil* b1) const;
};

/// A set of bilinear objects.
typedef std::set<LinBil*, CompareLinBil> LinBilSet;

/// Iterator of LinBil objects over a set.
typedef LinBilSet::iterator LinBilSetIter;
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
