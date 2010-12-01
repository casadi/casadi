/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#ifndef ACADO_INTERFACE_HPP
#define ACADO_INTERFACE_HPP

#include <casadi/fx/fx.hpp>

namespace CasADi{
  
  /** \brief  Inputs of an NLP Solver */
  enum ACADO_Input{
    ACADO_X_GUESS, // Initial guess for x [default: 0]
    ACADO_U_GUESS, // Initial guess for u [default: 0]
    ACADO_P_GUESS, // Initial guess for p [default: 0]
    ACADO_LBX, // Lower bound on x [default:  -infinity]
    ACADO_UBX, // Upper bound on x [default:  infinity]
    ACADO_LBX0, // Lower bound on x0 [default:  -infinity]
    ACADO_UBX0, // Upper bound on x0 [default:  infinity]
    ACADO_LBXF, // Lower bound on xf [default:  -infinity]
    ACADO_UBXF, // Upper bound on xf [default:  infinity]
    ACADO_LBU, // Lower bound on u [default:  -infinity]
    ACADO_UBU, // Upper bound on u [default:  infinity]
    ACADO_LBP, // Lower bound on p [default:  -infinity]
    ACADO_UBP, // Upper bound on p [default:  infinity]
    ACADO_LBC, // Lower bound on the path constraint function [default:  -infinity]
    ACADO_UBC, // Upper bound on the path constraint function [default:  infinity]
    ACADO_LBR, // Lower bound on the initial constraint function [default:  0]
    ACADO_UBR, // Upper bound on the initial constraint function [default:  0]
    ACADO_NUM_IN // Number of inputs
  };

  /** \brief  Outputs of an NLP Solver */
  enum ACADO_Output{
    ACADO_X_OPT,
    ACADO_U_OPT,
    ACADO_P_OPT,
    ACADO_COST,
    ACADO_NUM_OUT
  };

  /** \brief  Input arguments of an ACADO function */
  enum ACADO_FCN_Input{
    ACADO_FCN_T, 
    ACADO_FCN_XD, 
    ACADO_FCN_XA,
    ACADO_FCN_U,
    ACADO_FCN_P,
    ACADO_FCN_XDOT, 
    ACADO_FCN_NUM_IN
  };
  
  // Forward declaration
  class AcadoInternal;
 
  // Smart pointer class
  class AcadoInterface : public FX{
    public:

      /** \brief  Default constructor */
      AcadoInterface();

      /** \brief  Constructor taking a DAE rhs function, an objective function and a constraint function */
      explicit AcadoInterface(const FX& ffcn, const FX& mfcn, const FX& cfcn=FX(), const FX& rfcn=FX());

      /** \brief  Access functions and members of the node */
      AcadoInternal* operator->();

      /** \brief  Const access functions and members of the node */
      const AcadoInternal* operator->() const;
      
      /** \brief  Make sure that the pointer points towards a valid object */
      void assertNode() const;

};



} // namespace CasADi

#endif //ACADO_INTERFACE_HPP
