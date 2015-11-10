#include "extending_casadi.hpp"
#include <iostream>

namespace myproject {
  using namespace std;

  int ultimateQuestion() {
    return 42;
  }

  double casadi2scalar(const casadi::DM& x) {
    return x.toScalar();
  }

  casadi::DM scalar2casadi(double x) {
    return x;
  }

} // namespace myproject

