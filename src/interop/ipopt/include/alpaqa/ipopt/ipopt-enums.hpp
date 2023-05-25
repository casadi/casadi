#pragma once

#include <IpIpoptApplication.hpp>
#include <string_view>

inline std::string_view enum_name(Ipopt::ApplicationReturnStatus s) {
    using S = Ipopt::ApplicationReturnStatus;
    switch (s) {
        case S::Solve_Succeeded: return "Solve_Succeeded";
        case S::Solved_To_Acceptable_Level: return "Solved_To_Acceptable_Level";
        case S::Infeasible_Problem_Detected:
            return "Infeasible_Problem_Detected";
        case S::Search_Direction_Becomes_Too_Small:
            return "Search_Direction_Becomes_Too_Small";
        case S::Diverging_Iterates: return "Diverging_Iterates";
        case S::User_Requested_Stop: return "User_Requested_Stop";
        case S::Feasible_Point_Found: return "Feasible_Point_Found";
        case S::Maximum_Iterations_Exceeded:
            return "Maximum_Iterations_Exceeded";
        case S::Restoration_Failed: return "Restoration_Failed";
        case S::Error_In_Step_Computation: return "Error_In_Step_Computation";
        case S::Maximum_CpuTime_Exceeded: return "Maximum_CpuTime_Exceeded";
        case S::Maximum_WallTime_Exceeded: return "Maximum_WallTime_Exceeded";
        case S::Not_Enough_Degrees_Of_Freedom:
            return "Not_Enough_Degrees_Of_Freedom";
        case S::Invalid_Problem_Definition: return "Invalid_Problem_Definition";
        case S::Invalid_Option: return "Invalid_Option";
        case S::Invalid_Number_Detected: return "Invalid_Number_Detected";
        case S::Unrecoverable_Exception: return "Unrecoverable_Exception";
        case S::NonIpopt_Exception_Thrown: return "NonIpopt_Exception_Thrown";
        case S::Insufficient_Memory: return "Insufficient_Memory";
        case S::Internal_Error: return "Internal_Error";
        default: return "<unknown>";
    }
}
