// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file AMPLHessian.h
 * \brief Declare the AMPLHessian class for extracting Hessian of Lagrangian 
 * from AMPL.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURAMPLHESSIAN_H
#define MINOTAURAMPLHESSIAN_H

#include "asl.h"

#include "HessianOfLag.h"

namespace MINOTAUR_AMPL {
class AMPLInterface;
typedef AMPLInterface * AMPLInterfacePtr;

/**
 * \brief AMPLHessian class has methods to call Hessian of Lagrangian evaluation
 * routines of AMPL.  This class is not meant to calculate Hessian of native
 * Minotaur functions.
 */
class AMPLHessian : public Minotaur::HessianOfLag {
public:
  /**
   * \brief Construct AMPLHessian using an AMPL(ASL) interface.
   * \param[in] iface Pointer to AMPLInterface class from which we read the
   * instance.
   */
  AMPLHessian(AMPLInterfacePtr iface);

  /// Destroy.
  ~AMPLHessian();

  /// Return the number of non-zeros in the Hessian of Lagrangian.
  Minotaur::UInt getNumNz() const;

  /**
   * \brief Fill in row and column indices of nonzeros in Hessian of
   * Lagrangian.
   *
   * Given arrays irow and jcol, fill in the row and column index of each
   * non-zero in the Hessian of Lagrangian. The entries are ordered by
   * rows. 
   * \param [out] jcol Array of size returned by getNumNz().
   * \param [out] irow Array of size returned by getNumNz().
   */
  void fillRowColIndices(Minotaur::UInt *iRow, Minotaur::UInt *jCol);

  /** 
   * \brief Evaluate Hessian of Lagrangian at a given point x, row-wise.
   * \param[in] x The point at which Hessian of Lagrangian is to be evaluated.
   * \param[in] obj_mult The multiplier for the objective.
   * \param[in] con_mult Array of multipliers for constraints.
   * \param[out] values The array of size returned by getNumNz() into which
   * values are saved.
   * \param[out] error is set to a nonzero value if a problem is encountered, 
   * otherwise it is left undisturbed.
   */
  void fillRowColValues(const double *x, double obj_mult, const double
                        *con_mult, double *values, int *error);

  /// Ugly hack for maximization problems. 
  void negateObj();

private:
  /// Pointer to the AMPL's asl.
  ASL *myAsl_;

  /// No. of nonzeros in Hessian of Lagrangian. Initialized in the constructor.
  Minotaur::UInt nNz_;

  /// Size of tmp_ array.
  bool negObj_; 
};
typedef boost::shared_ptr<AMPLHessian> AMPLHessianPtr;
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

