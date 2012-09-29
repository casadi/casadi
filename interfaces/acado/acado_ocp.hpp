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

#include <symbolic/fx/fx.hpp>
#include <symbolic/fx/integrator.hpp>

namespace CasADi{
  
  /// Input arguments of an ACADO OCP solver [acadoIn]
  enum ACADO_Input{
    /// Initial guess for x (default: 0) [x_guess]
    ACADO_X_GUESS,
    /// Initial guess for u (default: 0) [u_guess]
    ACADO_U_GUESS,
    /// Initial guess for p (default: 0) [p_guess]
    ACADO_P_GUESS, 
    /// Lower bound on x (default:  -infinity) [lbx]
    ACADO_LBX,
    /// Upper bound on x (default:  infinity) [ubx]
    ACADO_UBX,
    /// Lower bound on x0 (default:  -infinity) [lbx0]
    ACADO_LBX0,
    /// Upper bound on x0 (default:  infinity) [ubx0]
    ACADO_UBX0,
    /// Lower bound on xf (default:  -infinity) [lbxf]
    ACADO_LBXF,
    /// Upper bound on xf (default:  infinity) [ubxf]
    ACADO_UBXF,
    /// Lower bound on u (default:  -infinity) [lbu]
    ACADO_LBU,
    /// Upper bound on u (default:  infinity) [ubu]
    ACADO_UBU,
    /// Lower bound on p (default:  -infinity) [lbp]
    ACADO_LBP,
    /// Upper bound on p (default:  infinity) [ubp]
    ACADO_UBP,
    /// Lower bound on the path constraint function (default:  -infinity) [lbc]
    ACADO_LBC,
    /// Upper bound on the path constraint function (default:  infinity) [ubc]
    ACADO_UBC,
    /// Lower bound on the initial constraint function (default:  0) [lbr]
    ACADO_LBR,
    /// Upper bound on the initial constraint function (default:  0) [ubr]
    ACADO_UBR,
    /// Number of inputs
    ACADO_NUM_IN
  };

  /// Output arguments of an ACADO OCP solver [acadoOut]
  enum ACADO_Output{
    /// Optimal states [x_opt]
    ACADO_X_OPT,
    /// Optimal control inputs [u_opt]
    ACADO_U_OPT,
    /// Optimal parameters [p_opt]
    ACADO_P_OPT,
    /// Optimal cost [cost]
    ACADO_COST,
    /// Number of outputs
    ACADO_NUM_OUT
  };

  /// Input arguments of an ACADO function [acadofcnIn]
  enum ACADO_FCN_Input{
    /// Time [t]
    ACADO_FCN_T,
    /// Differential state [xd] 
    ACADO_FCN_XD, 
    /// Algebraic state [xa]
    ACADO_FCN_XA,
    /// Control input [u]
    ACADO_FCN_U,
    /// Parameter [p]
    ACADO_FCN_P,
    /// Differential state derivative [xdot]
    ACADO_FCN_XDOT, 
    /// Number of inputs
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
