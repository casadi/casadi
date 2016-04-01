// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file AMPLJacobian.h
 * \brief Declare the AMPLJacobian class forreadin pcting Jacobian from AMPL.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURAMPLJACOBIAN_H
#define MINOTAURAMPLJACOBIAN_H

#include "asl.h"

#include "Jacobian.h"

namespace MINOTAUR_AMPL {

class AMPLInterface;
typedef AMPLInterface * AMPLInterfacePtr;

/**
 * \brief Jacobian class has methods to call Jacobian evaluation routines of AMPL.
 * This class is not meant to calculate Jacobian of native Minotaur
 * functions.
 */
class AMPLJacobian : public Minotaur::Jacobian {
public:
  /**
   * \brief Construct AMPLJacobian using an AMPL(ASL) interface.
   * \param[in] iface Pointer to AMPLInterface class from which we read the
   * instance.
   */
  AMPLJacobian(AMPLInterfacePtr iface);

  /// Destroy.
  ~AMPLJacobian();

  /**
   * \brief Fill in column and row indices of nonzeros in Jacobian.
   *
   * Given arrays irow and jcol, fill in the row and column index of each
   * non-zero in the jacobian. The entries are ordered by columns. This
   * is slightly inefficient because AMPL computes the gradients row-wise
   * and then translates them to column-wise (using goff).
   * e.g.
   * x1x3 + x1^2x4 = 2
   * x0 + x5 = 0
   * then, irow = [0 0 0 1 1]
   *       jcol = [1 3 4 0 5]
   *
   * \param [out] jcol Array of size returned by getNumNz().
   * \param [out] irow Array of size returned by getNumNz().
   */
  void fillColRowIndices(Minotaur::UInt *jcol, Minotaur::UInt *irow);

  /**
   * \brief Fill the indices of non-zeros in Jacobian, row-wise (or
   * constraint-wise).
   * \param [out] irow Array of size returned by getNumNz().
   * \param [out] jcol Array of size returned by getNumNz().
   */
  void fillRowColIndices(Minotaur::UInt *irow, Minotaur::UInt *jcol);

  /**
   * \brief Evaluate Jacobian at a given point x, column-wise.
   * Fill the values of Jacobian at a given point x into array 'values'. 
   * The ordering is same as that provided by
   * fillColRowIndices().
   * \param[in] x The point at which Jacobian is to be evaluated.
   * \param[out] values The array of size returned by getNumNz() into which
   * values are saved.
   * \param[out] error is set to a nonzero value if a problem is encountered,
   * otherwise it is left undisturbed.
   */
  void fillColRowValues(const double *x, 
                        double *values, int *error);

  /** 
   * \brief Evaluate Jacobian at a given point x, row-wise.
   * \param[in] x The point at which Jacobian is to be evaluated.
   * \param[out] values The array of size returned by getNumNz() into which
   * values are saved.
   * \param[out] error is set to a nonzero value if a problem is encountered, 
   * otherwise it is left undisturbed.
   */
  void fillRowColValues(const double *x, double *values, int *error);

  /// Return the number of non-zeros in the Jacobian.
  Minotaur::UInt getNumNz();

private:
  /// Pointer to the AMPL's asl.
  ASL *myAsl_;

  /// Number of nonzeros in Jacobian. Initialized in the constructor.
  Minotaur::UInt nNz_;

  /// Temporary array used to save Jacobian in fillRowColValues() function.
  double *tmp_;

  /// Size of tmp_ array.
  Minotaur::UInt tmpSize_;
};
typedef boost::shared_ptr<AMPLJacobian> AMPLJacobianPtr;
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

