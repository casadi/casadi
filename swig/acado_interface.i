%{
#include "acado_interface/acado_interface.hpp"
%}

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
    ACADO_LBC, // Lower bound on the constraint function [default:  -infinity]
    ACADO_UBC, // Upper bound on the constraint function [default:  infinity]
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

  class AcadoInterface : public FX{
    public:
    AcadoInterface(const FX& ffcn, const FX& mfcn, const FX& cfcn=FX(), const FX& rfcn=FX());
  };



} // namespace CasADi
