#include <alpaqa/alm.hpp>
#include <alpaqa/inner/directions/lbfgs.hpp>
#include <alpaqa/inner/panoc.hpp>

namespace alpaqa {
template class ALMSolver<PANOCSolver<LBFGS>>;
}