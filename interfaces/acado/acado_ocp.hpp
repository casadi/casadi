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

#ifndef ACADO_OCP_HPP
#define ACADO_OCP_HPP

#include <casadi/fx/fx.hpp>
#include <casadi/fx/integrator.hpp>

namespace CasADi{
  
  /// Inputs of an ACADO OCP solver
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

  /// Outputs of an ACADO OCP solver
  enum ACADO_Output{
    ACADO_X_OPT,
    ACADO_U_OPT,
    ACADO_P_OPT,
    ACADO_COST,
    ACADO_NUM_OUT
  };

  /// Input arguments of an ACADO function
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
  class AcadoOCPInternal;
 
  // Smart pointer class
  class AcadoOCP : public FX{
    public:

      /// Default constructor
      AcadoOCP();

      /// Constructor taking a DAE rhs function, an objective function and a constraint function -- for use with ACADO integrators
      explicit AcadoOCP(const FX& ffcn, const FX& mfcn, const FX& cfcn=FX(), const FX& rfcn=FX());

      /// Set a user-provided integrator
      void setIntegrators(const std::vector<Integrator>& integrators);
            
      /// Access functions and members of the node
      AcadoOCPInternal* operator->();

      /// Const access functions and members of the node
      const AcadoOCPInternal* operator->() const;
      
      /// Check if the node is pointing to the right type of object
      virtual bool checkNode() const;

};



} // namespace CasADi

#endif //ACADO_OCP_HPP
