#pragma once 

namespace alpaqa {

/// Which method to use to select the L-BFGS step size.
enum class LBFGSStepSize {
    BasedOnGradientStepSize = 0,
    BasedOnCurvature        = 1,
};

}