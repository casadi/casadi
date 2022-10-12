#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/inner/internal/panoc-stop-crit.hpp>
#include <alpaqa/inner/internal/solverstatus.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void register_enums(py::module_ &m) {

    py::enum_<alpaqa::LBFGSStepSize>(m, "LBFGSStepsize",
                                     "C++ documentation: :cpp:enum:`alpaqa::LBFGSStepSize`")
        .value("BasedOnExternalStepSize", alpaqa::LBFGSStepSize::BasedOnExternalStepSize)
        .value("BasedOnCurvature", alpaqa::LBFGSStepSize::BasedOnCurvature)
        .export_values();

    using SolverStatus = alpaqa::SolverStatus;
    py::enum_<SolverStatus>(m, "SolverStatus", py::arithmetic(),
                            "C++ documentation: :cpp:enum:`alpaqa::SolverStatus`")
        .value("Busy", SolverStatus::Busy, "In progress.")
        .value("Converged", SolverStatus::Converged, "Converged and reached given tolerance")
        .value("MaxTime", SolverStatus::MaxTime, "Maximum allowed execution time exceeded")
        .value("MaxIter", SolverStatus::MaxIter, "Maximum number of iterations exceeded")
        .value("NotFinite", SolverStatus::NotFinite, "Intermediate results were infinite or NaN")
        .value("NoProgress", SolverStatus::NoProgress, "No progress was made in the last iteration")
        .value("Interrupted", SolverStatus::Interrupted, "Solver was interrupted by the user")
        .export_values();

    py::enum_<alpaqa::PANOCStopCrit>(m, "PANOCStopCrit",
                                     "C++ documentation: :cpp:enum:`alpaqa::PANOCStopCrit`")
        .value("ApproxKKT", alpaqa::PANOCStopCrit::ApproxKKT)
        .value("ApproxKKT2", alpaqa::PANOCStopCrit::ApproxKKT2)
        .value("ProjGradNorm", alpaqa::PANOCStopCrit::ProjGradNorm)
        .value("ProjGradNorm2", alpaqa::PANOCStopCrit::ProjGradNorm2)
        .value("ProjGradUnitNorm", alpaqa::PANOCStopCrit::ProjGradUnitNorm)
        .value("ProjGradUnitNorm2", alpaqa::PANOCStopCrit::ProjGradUnitNorm2)
        .value("FPRNorm", alpaqa::PANOCStopCrit::FPRNorm)
        .value("FPRNorm2", alpaqa::PANOCStopCrit::FPRNorm2)
        .value("Ipopt", alpaqa::PANOCStopCrit::Ipopt)
        .value("LBFGSBpp", alpaqa::PANOCStopCrit::LBFGSBpp)
        .export_values();
}