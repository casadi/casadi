#pragma once

#include <alpaqa/accelerators/lbfgs.hpp>
#include <kwargs-to-struct.hpp>

PARAMS_TABLE_DECL(alpaqa::LBFGSParams<Conf>);
PARAMS_TABLE_DECL(alpaqa::CBFGSParams<Conf>);
