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

#ifndef SUNDIALS_INTERNAL_HPP
#define SUNDIALS_INTERNAL_HPP

#include "sundials_integrator.hpp"
#include "casadi/fx/integrator_internal.hpp"

namespace CasADi{
namespace Sundials{

class SundialsInternal : public IntegratorInternal{
public:
  /** \brief  Constructor */
  SundialsInternal(const FX& f, const FX& q);

  /** \brief  Destructor */
  virtual ~SundialsInternal()=0;
    
  /** \brief  Initialize */
  virtual void init();
};
  
} // namespace Sundials
} // namespace CasADi

#endif // SUNDIALS_INTERNAL_HPP
