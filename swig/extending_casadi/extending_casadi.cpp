#include "extending_casadi.hpp"
#include <iostream>

namespace myproject {

  int ultimateQuestion() {
    return 42;
  }

  double casadi2scalar(const casadi::DM& x) {
    return x.scalar();
  }

  casadi::DM scalar2casadi(double x) {
    return x;
  }

} // namespace myproject

