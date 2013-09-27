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

#ifndef GENERIC_INTEGRATOR_INTERNAL_HPP
#define GENERIC_INTEGRATOR_INTERNAL_HPP

#include "generic_integrator.hpp"
#include "fx_internal.hpp"

namespace CasADi{

  /** \brief Internal class for GenericIntegrator
      \author Joel Andersson 
      \date 2013
  */
  class GenericIntegratorInternal : public FXInternal{
  public:

    /** \brief  Constructor */
    GenericIntegratorInternal(const FX& f, const FX& g);

    /** \brief  Destructor */
    virtual ~GenericIntegratorInternal()=0;

    /** \brief  Clone */
    virtual GenericIntegratorInternal* clone() const=0;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);
  
    /** \brief  Create a new integrator */
    virtual GenericIntegratorInternal* create(const FX& f, const FX& g) const = 0;
  
    /** \brief  Initialize */
    virtual void init();
  
  };
  
} // namespace CasADi

#endif // GENERIC_INTEGRATOR_INTERNAL_HPP
