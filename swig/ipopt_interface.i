%{
#include "interfaces/ipopt/ipopt_solver.hpp"
#include <coin/IpReturnCodes.hpp>
%}

%include "interfaces/ipopt/ipopt_solver.hpp"



%inline %{
namespace CasADi {
  // Needed to expose in SWIG
  enum IPOPT_ApplicationReturnStatus {
    IPOPT_Solve_Succeeded = Ipopt::Solve_Succeeded,
    IPOPT_Solved_To_Acceptable_Level = Ipopt::Solved_To_Acceptable_Level,
    IPOPT_Infeasible_Problem_Detected = Ipopt::Infeasible_Problem_Detected,
    IPOPT_Search_Direction_Becomes_Too_Small = Ipopt::Search_Direction_Becomes_Too_Small,
    IPOPT_Diverging_Iterates = Ipopt::Diverging_Iterates,
    IPOPT_User_Requested_Stop = Ipopt::User_Requested_Stop,
    IPOPT_Maximum_Iterations_Exceeded = Ipopt::Maximum_Iterations_Exceeded,
    IPOPT_Restoration_Failed = Ipopt::Restoration_Failed,
    IPOPT_Error_In_Step_Computation = Ipopt::Error_In_Step_Computation,
    IPOPT_Not_Enough_Degrees_Of_Freedom = Ipopt::Not_Enough_Degrees_Of_Freedom,
    IPOPT_Invalid_Problem_Definition = Ipopt::Invalid_Problem_Definition,
    IPOPT_Invalid_Option = Ipopt::Invalid_Option,
    IPOPT_Invalid_Number_Detected = Ipopt::Invalid_Number_Detected,
    IPOPT_Unrecoverable_Exception = Ipopt::Unrecoverable_Exception,
    IPOPT_NonIpopt_Exception_Thrown = Ipopt::NonIpopt_Exception_Thrown,
    IPOPT_Insufficient_Memory = Ipopt::Insufficient_Memory
  };
}

%}
