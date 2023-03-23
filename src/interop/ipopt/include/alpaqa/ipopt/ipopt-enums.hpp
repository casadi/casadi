#pragma once

#include <IpIpoptApplication.hpp>
#include <string_view>

inline std::string_view enum_name(Ipopt::ApplicationReturnStatus s) {
    using enum Ipopt::ApplicationReturnStatus;
    switch (s) {
        case Solve_Succeeded: return "Solve_Succeeded";
        case Solved_To_Acceptable_Level: return "Solved_To_Acceptable_Level";
        case Infeasible_Problem_Detected: return "Infeasible_Problem_Detected";
        case Search_Direction_Becomes_Too_Small:
            return "Search_Direction_Becomes_Too_Small";
        case Diverging_Iterates: return "Diverging_Iterates";
        case User_Requested_Stop: return "User_Requested_Stop";
        case Feasible_Point_Found: return "Feasible_Point_Found";
        case Maximum_Iterations_Exceeded: return "Maximum_Iterations_Exceeded";
        case Restoration_Failed: return "Restoration_Failed";
        case Error_In_Step_Computation: return "Error_In_Step_Computation";
        case Maximum_CpuTime_Exceeded: return "Maximum_CpuTime_Exceeded";
        case Maximum_WallTime_Exceeded: return "Maximum_WallTime_Exceeded";
        case Not_Enough_Degrees_Of_Freedom:
            return "Not_Enough_Degrees_Of_Freedom";
        case Invalid_Problem_Definition: return "Invalid_Problem_Definition";
        case Invalid_Option: return "Invalid_Option";
        case Invalid_Number_Detected: return "Invalid_Number_Detected";
        case Unrecoverable_Exception: return "Unrecoverable_Exception";
        case NonIpopt_Exception_Thrown: return "NonIpopt_Exception_Thrown";
        case Insufficient_Memory: return "Insufficient_Memory";
        case Internal_Error: return "Internal_Error";
        default: return "<unknown>";
    }
}
