#ifndef EXTENDING_CASADI_HPP
#define EXTENDING_CASADI_HPP

#include <casadi/casadi.hpp>

namespace myproject {
  
  /** \brief Answer to the ultimate question to life, the universe and everything
      You know the answer.
  */
  int ultimateQuestion();

  /** \brief Convert a CasADi DM to a scalar
      Argument must be a 1-by-1 matrix.
  */
  double casadi2scalar(const casadi::DM& x);

  /** \brief Convert a scalar to a CasADi DM
      Argument must be a 1-by-1 matrix.
  */
  casadi::DM scalar2casadi(double x);
  
} // namespace myproject

#endif // EXTENDING_CASADI_HPP
