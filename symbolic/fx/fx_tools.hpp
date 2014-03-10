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

#ifndef FX_TOOLS_HPP
#define FX_TOOLS_HPP

#include "fx.hpp"
#include "mx_function.hpp"

namespace CasADi{

    /** \brief make integration start and end time a parameter
    Applies the conversion  t = t0 + (tf-t0)*tau to the supplied dae.
    with tau dimensionless time.
    The Input/OuputScheme of the result is the same as the scheme of the dae, except for input(DAE_P), which is extended by t0 and tf at the top.
    */
    FX parameterizeTime(FX dae);
    
    /** \brief adapts an output function such that start and end time are parameters
    Applies the conversion  t = t0 + (tf-t0)*tau to the supplied dae.
    with tau dimensionless time.
    The InputScheme of the result is the same as the scheme of the dae, except for input(DAE_P), which is extended by t0 and tf at the top.
    */
    FX parameterizeTimeOutput(FX outputfcn);
                            
} // namespace CasADi


#endif // FX_TOOLS_HPP
