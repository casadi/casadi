#pragma once

#include <alpaqa/inner/directions/panoc/structured-newton.hpp>
#include <kwargs-to-struct.hpp>

PARAMS_TABLE_DECL(alpaqa::StructuredNewtonRegularizationParams<Conf>);
PARAMS_TABLE_DECL(alpaqa::StructuredNewtonDirectionParams<Conf>);
