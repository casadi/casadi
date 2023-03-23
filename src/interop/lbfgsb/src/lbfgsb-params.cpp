#include <alpaqa/implementation/params/params.tpp>
#include <alpaqa/lbfgsb/lbfgsb-adapter.hpp>

namespace alpaqa::params {

PARAMS_TABLE(lbfgsb::LBFGSBSolver::Params,   //
             PARAMS_MEMBER(memory),          //
             PARAMS_MEMBER(max_iter),        //
             PARAMS_MEMBER(max_time),        //
             PARAMS_MEMBER(stop_crit),       //
             PARAMS_MEMBER(print),           //
             PARAMS_MEMBER(print_interval),  //
             PARAMS_MEMBER(print_precision), //
);

template void LBFGSB_ADAPTER_EXPORT
set_params(lbfgsb::LBFGSBSolver::Params &, std::string_view,
           std::span<const std::string_view>, std::optional<std::span<bool>>);

} // namespace alpaqa::params