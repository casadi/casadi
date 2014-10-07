/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_SMITH_INDEF_DLE_INTERNAL_HPP
#define CASADI_SMITH_INDEF_DLE_INTERNAL_HPP

#include "../core/function/lr_dle_internal.hpp"
#include <casadi/solvers/casadi_lrdlesolver_smith_export.h>

/** \defgroup LrDleSolversmith

   For the low rank case, we have:
   
   \verbatim
       P0 = CVC^T
       P1 = ACVC^TA^T + CVC^T
       P2 = AACVC^TAA^TAA^T + ACVC^TA^T + CVC^T
       .....
   \endverbatim


   In other words, in each iteration, we perform a low-rank update
   to the initial value of P.

         P_k = P_{k-1} + D_k V D_k

   with

   D = [ C AC AAC AAAC ... ]
   
   There is no need to actually store D:
   
   \verbatim
     C_0 = C
     Y_0 = 0
     k = 0
     
     while || C_k V C_k^T || < epsilon
       Y_{k+1} = Y_k + H^T C_k V C_k^T H
       C_{k+1} = A_k C_k
       k += 1
     end
     
     Y = Y_k
   \endverbatim
*/

/** \defgroup plugin_LrDleSolver_smith
 Solving the Low-rank Discrete Lyapunov Equations with Smith iterations

   @copdyoc DleSolversmith
   @copdyoc LrDleSolversmith


   Implementation details:
     - We avoid ever holding P in memory as it might be large
       norm_inf_mul_tt was used to obtain a stopping criteria

     - We avoid memory allocation in evaluate.
       All sparsity pattern calculations have been done at init
    
*/
/** \pluginsection{LrDleSolver,smith} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{LrDleSolver,smith}

   @copydoc LR_DLE_doc
   @copydoc plugin_LrDleSolver_smith

       \author Joris Gillis
      \date 2014

  */
  class CASADI_LRDLESOLVER_SMITH_EXPORT SmithLrDleInternal : public LrDleInternal {
  public:
    /** \brief  Constructor
     * \param st \structargument{LrDle}
     * \param Hs Column-sizes of H_i
     */
    SmithLrDleInternal(const LrDleStructure& st, const std::vector<int> &Hs,
      int nfwd=0, int nadj=0);

    /** \brief  Destructor */
    virtual ~SmithLrDleInternal();

    /** \brief  Clone */
    virtual SmithLrDleInternal* clone() const;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new solver */
    virtual SmithLrDleInternal* create(const LrDleStructure& st, const std::vector<int> &Hs) const {
        return new SmithLrDleInternal(st, Hs);}

    /** \brief  Create a new DLE Solver */
    static LrDleInternal* creator(const LrDleStructure& st, const std::vector<int> &Hs)
    { return new SmithLrDleInternal(st, Hs);}

    /** \brief  Print solver statistics */
    virtual void printStats(std::ostream &stream) const {}

    /** \brief  evaluate */
    virtual void evaluate();

    /** \brief  Initialize */
    virtual void init();

    /** \brief Generate a function that calculates \a nfwd forward derivatives
     and \a nadj adjoint derivatives
    */
    virtual Function getDerivative(int nfwd, int nadj);

    /// A documentation string
    static const std::string meta_doc;

  private:

    /// State space dimension
    int n_;

    /// If true, each iteration will be printed
    bool print_iteration_;

    /// Print iteration header
    void printIteration(std::ostream &stream);

    /// Print iteration
    void printIteration(std::ostream &stream, int iter, double norm_inf);

    // Forward seed DLE_V symmetrized
    DMatrix Vdsymm_;

    // Input DLE_V symmetrized
    DMatrix Vsymm_;

    /// Intermediate variables
    std::vector<DMatrix> D_;   // Contains AAAAA...C
    std::vector<DMatrix> VD_;  // V D^T
    std::vector<DMatrix> DTH_; // D^T H
    std::vector< std::vector<DMatrix> > DTHv_;
    std::vector< std::vector<DMatrix> > VDTHv_;
    std::vector<DMatrix> VDTH_; // V D^T H

    /// Tape for AD
    std::vector<DMatrix> D_ad_;
    std::vector<DMatrix> VD_ad_;
    std::vector<DMatrix> DTH_ad_;
    std::vector< std::vector<DMatrix> > DTHv_ad_;
    std::vector< std::vector<DMatrix> > VDTHv_ad_;
    std::vector<DMatrix> VDTH_ad_;

    // Indices to effectively slice H_j from H
    std::vector< std::vector<int> > DTHvi_;
    std::vector< std::vector<int> > VDTHvi_;

    DMatrix e_;
    DMatrix em_;

    DMatrix CV_;

    // Work vectors for norm_inf_mul_tt
    std::vector<double> Dwork_norm_;
    std::vector<int> Iwork_norm_;

    /// Maximum number of iterations for the algorithm
    int max_iter_;

    /// Tolerance for satisfying the Lyapunov equation
    double tol_;

    /// Number of forward directions
    int nfwd_;

    /// Number of adjoint directions
    int nadj_;

    /// Prepare intermediate variables and AD tape
    void smith_prepare();

  };

} // namespace casadi
/// \endcond
#endif // CASADI_SMITH_INDEF_DLE_INTERNAL_HPP
