// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file Operations.h
 * \brief Declare some commonly used functions.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef OPERATIONS_H
#define OPERATIONS_H

#include "LinearFunction.h"

namespace Minotaur {
  bool IsInt(double v, double tol=1e-12);

  double InnerProduct(const VariableGroup &g1, const VariableGroup &g2);

  double InnerProduct(const std::vector<double> &x, const VariableGroup &g2);

  /// Inner product of two arrays, x and a, each of size n.
  double InnerProduct(const double *x, const double *a, int n);

  void symMatDotV(UInt nz, const double *mat, const UInt *irow, const UInt *jcol,
                  const double *v, double *prod);

  void BoundsOnDiv(double l0, double u0, double l1, double u1, 
                   double &lb, double &ub);

  /// Get bounds on y, where \f$ y = x_0x_1 \f$.
  void BoundsOnProduct(ConstVariablePtr x0, ConstVariablePtr x1,
                       double &lb, double &ub);

  /// Get bounds on y, where \f$ y = x_0x_1 \f$.
  void BoundsOnProduct(double l0, double u0, double l1, double u1, double &lb,
                       double &ub);

  void BoundsOnRecip(double l0, double u0, double &lb, double &ub);

  /// Get lower and upper bounds on a square of variable.
  void BoundsOnSquare(ConstVariablePtr x1, double &lb, double &ub);

  /**
   * Get lower and upper bounds on a square of variable whose bounds are
   * provided.
   */
  void BoundsOnSquare(const double x_lb, const double y_lb, double &lb, 
                      double &ub);

  /**
   * Get gcd(d1, d2). Returns d2 if d1 is zero. Returns d1 if d2 is zero.
   * Otherwise, returns largest x > 0 s.t. d1/x and d2/x are integers.
   * e.g. Gcd(10.5, 7.0) = 3.5. If d1 is zero, returns abs(d2), if d2 is zero,
   * returns abs(d1).
   */
  double Gcd(double d1, double d2, const double &etol = 1e-9);

  /**
   * Get distance(a,b,n). Returns the Euclidian distance between point a and point b which 
   * in n-dimensional space.
   */
  double getDistance(const double* Pointa, const double* Pointb, UInt n);

  /**
   * Get minArray. Returns the minimum element of a dynamic array.
   */
  double minArray(const double* A, UInt n);

  /// To display the initial point.
  void displayArray(const double* point, UInt n, std::ostream &out);

  /// Sort a vector of variables according to values of x.
  void sort(VarVector &vvec, double *x, bool ascend=true);
  void sortRec(VarVector &vvec, double *x, int left, int right, int pivotind);

  /// Convert a string to lower case.
  void toLowerCase(std::string &str);

  /// Convert a double representing seconds into wall clock time HH:MM:SS:mm
  std::string toClockTime(double t);
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
