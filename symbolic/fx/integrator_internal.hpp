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

#ifndef INTEGRATOR_INTERNAL_HPP
#define INTEGRATOR_INTERNAL_HPP

#include "integrator.hpp"
#include "generic_integrator_internal.hpp"

namespace CasADi{

  /** \brief Internal storage for integrator related data

      @copydoc DAE_doc
      \author Joel Andersson 
      \date 2010
  */
  class IntegratorInternal : public GenericIntegratorInternal{
  public:
    /** \brief  Constructor */
    IntegratorInternal(const FX& f, const FX& g);

    /** \brief  Destructor */
    virtual ~IntegratorInternal()=0;

    /** \brief  Clone */
    virtual IntegratorInternal* clone() const=0;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);
  
    /** \brief  Create a new integrator */
    virtual IntegratorInternal* create(const FX& f, const FX& g) const = 0;
  
    /** \brief  Print solver statistics */
    virtual void printStats(std::ostream &stream) const{}

    /** \brief  Reset the forward problem and bring the time back to t0 */
    virtual void reset() = 0;

    /** \brief  Reset the backward problem and take time to tf */
    virtual void resetB() = 0;

    /** \brief  Integrate forward until a specified time point */
    virtual void integrate(double t_out) = 0;

    /** \brief  Integrate backward until a specified time point */
    virtual void integrateB(double t_out) = 0;

    /** \brief  evaluate */
    virtual void evaluate();

    /** \brief  Initialize */
    virtual void init();

    /** \brief  Propagate the sparsity pattern through a set of directional derivatives forward or backward */
    virtual void spEvaluate(bool fwd);

    /// Is the class able to propate seeds through the algorithm?
    virtual bool spCanEvaluate(bool fwd){ return true;}
    
    /// Generate a function that calculates nfwd forward derivatives and nadj adjoint derivatives
    virtual FX getDerivative(int nfwd, int nadj);

    /** \brief Calculate the jacobian of output oind with respect to input iind */
    virtual FX getJacobian(int iind, int oind, bool compact, bool symmetric);

    /// Generate a augmented DAE system with nfwd forward sensitivities and nadj adjoint sensitivities
    virtual std::pair<FX,FX> getAugmented(int nfwd, int nadj, std::vector<int>& xf_offset, std::vector<int>& qf_offset, std::vector<int>& rxf_offset, std::vector<int>& rqf_offset);
  
    /// Generate a augmented DAE system with nfwd forward sensitivities and nadj adjoint sensitivities (generic)
    template<class Mat,class XFunc>
    std::pair<FX,FX> getAugmentedGen(int nfwd, int nadj, std::vector<int>& xf_offset, std::vector<int>& qf_offset, std::vector<int>& rxf_offset, std::vector<int>& rqf_offset);
  
    /// Integration horizon
    double t0_, tf_;
  
    /// ODE/DAE forward integration function
    FX f_;
  
    /// ODE/DAE backward integration function, if any
    FX g_;
    
    /// Algebraic variable
    DMatrix z_, rz_;
  };
  
} // namespace CasADi

#endif // INTEGRATOR_INTERNAL_HPP
