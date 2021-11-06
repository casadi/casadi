#include <alpaqa/alm.hpp>
#include <alpaqa/inner/directions/lbfgs.hpp>
#include <alpaqa/inner/panoc.hpp>

namespace pa {
template class ALMSolver<PANOCSolver<LBFGS>>;
}