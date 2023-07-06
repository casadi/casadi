#ifndef SLEQP_FUNC_H
#define SLEQP_FUNC_H

#include <sleqp.h>

namespace casadi {
  struct SLEQPMemory;

  void casadi_sleqp_func_create(SleqpFunc** star,
                                int num_vars,
                                int num_cons,
                                SLEQPMemory* m);

} // namespace casadi

#endif /* SLEQP_FUNC_H */
