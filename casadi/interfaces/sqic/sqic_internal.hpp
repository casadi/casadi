/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef SQIC_INTERNAL_HPP
#define SQIC_INTERNAL_HPP

#include "casadi/core/function/qp_solver_internal.hpp"
#include "sqic_solver.hpp"

/// \cond INTERNAL
namespace casadi {

/** \brief Internal class for SQICSolver
 *
    @copydoc QPSolver_doc
 * */
class CASADI_SQIC_INTERFACE_EXPORT SQICInternal : public QPSolverInternal {
  friend class SQICSolver;

public:
  /** \brief  Constructor */
  explicit SQICInternal();

  /** \brief  Clone */
  virtual SQICInternal* clone() const;

  /** \brief  Create a new Solver */
  explicit SQICInternal(const std::vector<Sparsity>& st);

  /** \brief  Destructor */
  virtual ~SQICInternal();

  /** \brief  Initialize */
  virtual void init();

  /** \brief Generate native code for debugging */
  virtual void generateNativeCode(std::ostream& file) const;

  virtual void evaluate();

    /// Throw error
    static void sqic_error(const std::string& module, int flag);

    /// Calculate the error message map
    static std::map<int, std::string> calc_flagmap();

    /// Error message map
    static std::map<int, std::string> flagmap;
  protected:

    /// Flag: is already initialized
    bool is_init_;

    /// Storage space for sqic \p bl variable
    std::vector<double> bl_;
    /// Storage space for sqic \p bu variable
    std::vector<double> bu_;
    /// Storage space for sqic \p x variable
    std::vector<double> x_;
    /// Storage space for sqic \p locA variable
    std::vector<int> locA_;
    /// Storage space for sqic \p indA variable
    std::vector<int> indA_;
    /// Storage space for sqic \p hs variable
    std::vector<int> hs_;
    /// Storage space for sqic \p hEtype variable
    std::vector<int> hEtype_;
    /// Storage space for sqic \p indH variable
    std::vector<int> indH_;
    /// Storage space for sqic \p locH variable
    std::vector<int> locH_;
    /// Storage space for sqic \p rc variable
    std::vector<double> rc_;
    /// Storage space for sqic \p rc variable
    std::vector<double> pi_;

    /// Helper function to bring A into correct format
    Function formatA_;

    /// SQIC inf
    double inf_;


};


} // namespace casadi

/// \endcond
#endif //SQIC_INTERNAL_HPP

