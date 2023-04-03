//
//    MIT No Attribution
//
//    Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

typedef struct {
  // Name of instance
  const char* instance_name;
  // Value for all variables
  double v[SZ_MEM];
  // Value for derivatives
  double d[SZ_MEM];
  // Has any variable been set since last call
  int up_to_date;
  // Buffers for evaluation generated code
  double t;
  double p[N_P];
  double x[N_X];
  double xdot[N_X];
  double u[N_U];
  double y[N_Y];
  // Buffers for derivative calculations
  double dp[N_P];
  double dx[N_X];
  double dxdot[N_X];
  double du[N_U];
  double dy[N_Y];
  // Work vectors for evaluation
  const double* arg[SZ_ARG];
  double* res[SZ_RES];
  casadi_int iw[SZ_IW];
  double w[SZ_W];
} casadi_fmi_memory;

int evaluate(casadi_fmi_memory* m) {
  // Local variables
  size_t i;
  int mem, flag;
  // Copy states, inputs and parameters to input buffers
  for (i = 0; i < N_X; ++i) m->x[i] = m->v[x_vr[i]];
  for (i = 0; i < N_P; ++i) m->p[i] = m->v[p_vr[i]];
  for (i = 0; i < N_U; ++i) m->u[i] = m->v[u_vr[i]];

  // Map inputs to evaluation buffer
  i = 0;
  m->arg[i++] = &m->t;
  m->arg[i++] = m->x;
  m->arg[i++] = m->p;
  m->arg[i++] = m->u;

  // Map outputs to evaluation buffer
  i = 0;
  m->res[i++] = m->xdot;
  m->res[i++] = m->y;

  // Evaluate
  mem = MODELNAME_checkout();
  flag = MODELNAME(m->arg, m->res, m->iw, m->w, mem);
  MODELNAME_release(mem);

  // Copy from output buffers
  for (i = 0; i < N_X; ++i) m->v[xdot_vr[i]] = m->xdot[i];
  for (i = 0; i < N_Y; ++i) m->v[y_vr[i]] = m->y[i];

  return flag;
}

int evaluate_forward(casadi_fmi_memory* m) {
  // Local variables
  size_t i;
  int mem, flag;
  // Copy seeds for states, inputs and parameters to input buffers
  for (i = 0; i < N_X; ++i) m->dx[i] = m->d[x_vr[i]];
  for (i = 0; i < N_P; ++i) m->dp[i] = m->d[p_vr[i]];
  for (i = 0; i < N_U; ++i) m->du[i] = m->d[u_vr[i]];

  // Map nondifferentiated inputs to evaluation buffer
  i = 0;
  m->arg[i++] = 0;  // t
  m->arg[i++] = m->x;
  m->arg[i++] = m->p;
  m->arg[i++] = m->u;

  // Map nondifferentiated outputs to evaluation buffer
  m->arg[i++] = m->xdot;
  m->arg[i++] = m->y;

  // Map forward seeds
  m->arg[i++] = 0;
  m->arg[i++] = m->dx;
  m->arg[i++] = m->dp;
  m->arg[i++] = m->du;

  // Map forward sensitivities
  i = 0;
  m->res[i++] = m->dxdot;
  m->res[i++] = m->dy;

  // Evaluate
  mem = fwd1_MODELNAME_checkout();
  flag = fwd1_MODELNAME(m->arg, m->res, m->iw, m->w, mem);
  fwd1_MODELNAME_release(mem);

  // Copy from output buffers
  for (i = 0; i < N_X; ++i) m->d[xdot_vr[i]] = m->dxdot[i];
  for (i = 0; i < N_Y; ++i) m->d[y_vr[i]] = m->dy[i];

  // Return evaluation flag
  return flag;
}

int evaluate_adjoint(casadi_fmi_memory* m) {
  // Local variables
  size_t i;
  int mem, flag;
  // Copy seeds for output buffers
  for (i = 0; i < N_X; ++i) m->dxdot[i] = m->d[xdot_vr[i]];
  for (i = 0; i < N_Y; ++i) m->dy[i] = m->d[y_vr[i]];

  // Clear sensitivities
  for (i = 0; i < N_X; ++i) m->dx[i] = 0;
  for (i = 0; i < N_P; ++i) m->dp[i] = 0;
  for (i = 0; i < N_U; ++i) m->du[i] = 0;

  // Map nondifferentiated inputs to evaluation buffer
  i = 0;
  m->arg[i++] = 0;  // t
  m->arg[i++] = m->x;
  m->arg[i++] = m->p;
  m->arg[i++] = m->u;

  // Map nondifferentiated outputs to evaluation buffer
  m->arg[i++] = m->xdot;
  m->arg[i++] = m->y;

  // Map adjoint seeds
  m->arg[i++] = m->dxdot;
  m->arg[i++] = m->dy;

  // Map adjoint sensitivities
  i = 0;
  m->res[i++] = 0;  // t
  m->res[i++] = m->dx;
  m->res[i++] = m->dp;
  m->res[i++] = m->du;

  // Evaluate
  mem = adj1_MODELNAME_checkout();
  flag = adj1_MODELNAME(m->arg, m->res, m->iw, m->w, mem);
  adj1_MODELNAME_release(mem);

  // Copy from input buffers
  for (i = 0; i < N_X; ++i) m->d[x_vr[i]] = m->dx[i];
  for (i = 0; i < N_P; ++i) m->d[p_vr[i]] = m->dp[i];
  for (i = 0; i < N_U; ++i) m->d[u_vr[i]] = m->du[i];

  // Return evaluation flag
  return flag;
}

FMI3_Export fmi3Status fmi3Reset(fmi3Instance instance) {
  // Local variables
  size_t i;
  casadi_fmi_memory* m;
  // Cast to memory struct
  m = (casadi_fmi_memory*)instance;
  // Initialize variables
  for (i = 0; i < SZ_MEM; ++i) m->v[i] = start[i];
  // Reset derivative buffer
  for (i = 0; i < SZ_MEM; ++i) m->d[i] = 0;
  // Not evaluated
  m->up_to_date = 0;
  // Always successful return
  return fmi3OK;
}

FMI3_Export fmi3Instance fmi3InstantiateModelExchange(
    fmi3String instanceName,
    fmi3String instantiationToken,
    fmi3String resourcePath,
    fmi3Boolean visible,
    fmi3Boolean loggingOn,
    fmi3InstanceEnvironment instanceEnvironment,
    fmi3LogMessageCallback logMessage) {
  // Local variables
  casadi_fmi_memory* m;
  // Unused variables
  (void)resourcePath;  // unused
  (void)visible;  // unused
  (void)loggingOn;  // unused
  (void)instanceEnvironment;  // unused
  (void)logMessage;  // unused
  // Allocate memory structure
  m = (casadi_fmi_memory*)malloc(sizeof(casadi_fmi_memory));
  // If allocation was successful
  if (m) {
    // Increase counter for codegen
    MODELNAME_incref();
    // Copy meta data
    m->instance_name = instanceName;
    // Call reset function (return flag does not need to be checked)
    (void)fmi3Reset(m);
  }
  // Return pointer to instance, or null
  return m;
}

FMI3_Export void fmi3FreeInstance(fmi3Instance instance) {
  if (instance) {
    // Free memory structure
    free(instance);
    // Decrease counter for codegen
    MODELNAME_decref();
  }
}

FMI3_Export fmi3Status fmi3SetFloat64(
    fmi3Instance instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    const fmi3Float64 values[],
    size_t nValues) {
  // Local variables
  casadi_fmi_memory* m;
  size_t i, j, var_off, var_sz, val_ind;
  fmi3ValueReference vr;
  // Cast to memory struct
  m = (casadi_fmi_memory*)instance;
  // Not evaluated
  m->up_to_date = 0;
  // Position in values vector
  val_ind = 0;
  // Loop over provided variables
  for (i = 0; i < nValueReferences; ++i) {
    // Get variable
    vr = valueReferences[i];
    // Get offset, size of variable
    var_off = var_offset[vr];
    var_sz = var_offset[vr + 1] - var_off;
    // Copy elements
    for (j = 0; j < var_sz; ++j) m->v[var_off + j] = values[val_ind++];
  }
  // Consistency check
  if (val_ind != nValues) return fmi3Fatal;
  // Successful return
  return fmi3OK;
}

FMI3_Export fmi3Status fmi3GetFloat64(
    fmi3Instance instance,
    const fmi3ValueReference valueReferences[],
    size_t nValueReferences,
    fmi3Float64 values[],
    size_t nValues) {
  // Local variables
  casadi_fmi_memory* m;
  size_t i, j, var_off, var_sz, val_ind;
  fmi3ValueReference vr;
  // Cast to memory struct
  m = (casadi_fmi_memory*)instance;
  // Evaluate, if necessary
  if (!m->up_to_date) {
    // Evaluate state derivatives and outputs
    if (evaluate(m)) return fmi3Error;
    // Now evaluated
    m->up_to_date = 1;
  }
  // Position in values vector
  val_ind = 0;
  // Loop over provided variables
  for (i = 0; i < nValueReferences; ++i) {
    // Get variable
    vr = valueReferences[i];
    // Get offset, size of variable
    var_off = var_offset[vr];
    var_sz = var_offset[vr + 1] - var_off;
    // Copy elements
    for (j = 0; j < var_sz; ++j) values[val_ind++] = m->v[var_off + j];
  }
  // Consistency check
  if (val_ind != nValues) return fmi3Fatal;
  // Successful return
  return fmi3OK;  
}

FMI3_Export fmi3Status fmi3EnterInitializationMode(
    fmi3Instance instance,
    fmi3Boolean toleranceDefined,
    fmi3Float64 tolerance,
    fmi3Float64 startTime,
    fmi3Boolean stopTimeDefined,
    fmi3Float64 stopTime) {
  // Local variables
  casadi_fmi_memory* m;
  // Unused variables
  (void)toleranceDefined;  // unused
  (void)tolerance;  // unused
  (void)stopTimeDefined;  // unused
  (void)stopTime;  // unused
  // Cast to memory struct
  m = (casadi_fmi_memory*)instance;
  // Initialize time
  m->t = startTime;
  // Always successful return
  return fmi3OK;
}

FMI3_Export fmi3Status fmi3ExitInitializationMode(fmi3Instance instance) {
  // Unused variables
  (void)instance;  // unused
  // Always successful return
  return fmi3OK;
}

FMI3_Export fmi3Status fmi3EnterContinuousTimeMode(fmi3Instance instance) {
  // Unused variables
  (void)instance;  // unused
  // Always successful return
  return fmi3OK;
}

FMI3_Export fmi3Status fmi3SetTime(fmi3Instance instance, fmi3Float64 time) {
  // Local variables
  casadi_fmi_memory* m;
  // Cast to memory struct
  m = (casadi_fmi_memory*)instance;
  // Initialize time
  m->t = time;
  // Always successful return
  return fmi3OK;
}

FMI3_Export fmi3Status fmi3SetContinuousStates(
    fmi3Instance instance,
    const fmi3Float64 continuousStates[],
    size_t nContinuousStates) {
  return fmi3SetFloat64(instance, x_vr, N_X, continuousStates, nContinuousStates);
}

FMI3_Export fmi3Status fmi3GetContinuousStates(
    fmi3Instance instance,
    fmi3Float64 continuousStates[],
    size_t nContinuousStates) {
  return fmi3GetFloat64(instance, x_vr, N_X, continuousStates, nContinuousStates);
}

fmi3Status fmi3GetDirectionalDerivative(
    fmi3Instance instance,
    const fmi3ValueReference unknowns[],
    size_t nUnknowns,
    const fmi3ValueReference knowns[],
    size_t nKnowns,
    const fmi3Float64 seed[],
    size_t nSeed,
    fmi3Float64 sensitivity[],
    size_t nSensitivity) {
  // Local variables
  casadi_fmi_memory* m;
  size_t i, j, var_off, var_sz, val_ind;
  fmi3ValueReference vr;
  int flag;
  // Cast to memory struct
  m = (casadi_fmi_memory*)instance;
  // Evaluate, if necessary
  if (!m->up_to_date) {
    // Evaluate state derivatives and outputs
    if (evaluate(m)) return fmi3Error;
    // Now evaluated
    m->up_to_date = 1;
  }
  // Pass derivative seeds
  val_ind = 0;
  for (i = 0; i < nKnowns; ++i) {
    // Get variable
    vr = knowns[i];
    // Get offset, size of variable
    var_off = var_offset[vr];
    var_sz = var_offset[vr + 1] - var_off;
    // Copy elements
    for (j = 0; j < var_sz; ++j) m->d[var_off + j] = seed[val_ind++];
  }
  // Consistency check
  if (val_ind != nSeed) return fmi3Fatal;
  // Evaluate forward directional derivatives
  flag = evaluate_forward(m);
  // Collect sensitivities
  val_ind = 0;
  for (i = 0; i < nUnknowns; ++i) {
    // Get variable
    vr = unknowns[i];
    // Get offset, size of variable
    var_off = var_offset[vr];
    var_sz = var_offset[vr + 1] - var_off;
    // Copy elements, clear d
    for (j = 0; j < var_sz; ++j) {
      sensitivity[val_ind++] = m->d[var_off + j];
      m->d[var_off + j] = 0;
    }
  }
  // Consistency check
  if (val_ind != nSensitivity) return fmi3Fatal;
  // Clear derivative seeds
  val_ind = 0;
  for (i = 0; i < nKnowns; ++i) {
    // Get variable
    vr = knowns[i];
    // Get offset, size of variable
    var_off = var_offset[vr];
    var_sz = var_offset[vr + 1] - var_off;
    // Clear elements
    for (j = 0; j < var_sz; ++j) m->d[var_off + j] = 0;
  }
  // Check for evaluation error
  if (flag) return fmi3Error;
  // Successful return
  return fmi3OK;
}

fmi3Status fmi3GetAdjointDerivative(
    fmi3Instance instance,
    const fmi3ValueReference unknowns[],
    size_t nUnknowns,
    const fmi3ValueReference knowns[],
    size_t nKnowns,
    const fmi3Float64 seed[],
    size_t nSeed,
    fmi3Float64 sensitivity[],
    size_t nSensitivity) {
  // Local variables
  casadi_fmi_memory* m;
  size_t i, j, var_off, var_sz, val_ind;
  fmi3ValueReference vr;
  int flag;
  // Cast to memory struct
  m = (casadi_fmi_memory*)instance;
  // Evaluate, if necessary
  if (!m->up_to_date) {
    // Evaluate state derivatives and outputs
    if (evaluate(m)) return fmi3Error;
    // Now evaluated
    m->up_to_date = 1;
  }
  // Pass derivative seeds
  val_ind = 0;
  for (i = 0; i < nUnknowns; ++i) {
    // Get variable
    vr = unknowns[i];
    // Get offset, size of variable
    var_off = var_offset[vr];
    var_sz = var_offset[vr + 1] - var_off;
    // Copy elements
    for (j = 0; j < var_sz; ++j) m->d[var_off + j] = seed[val_ind++];
  }
  // Consistency check
  if (val_ind != nSeed) return fmi3Fatal;
  // Evaluate adjoint directional derivatives
  flag = evaluate_adjoint(m);
  // Collect sensitivities
  val_ind = 0;
  for (i = 0; i < nKnowns; ++i) {
    // Get variable
    vr = knowns[i];
    // Get offset, size of variable
    var_off = var_offset[vr];
    var_sz = var_offset[vr + 1] - var_off;
    // Copy elements, clear d
    for (j = 0; j < var_sz; ++j) {
      sensitivity[val_ind++] = m->d[var_off + j];
      m->d[var_off + j] = 0;
    }
  }
  // Consistency check
  if (val_ind != nSensitivity) return fmi3Fatal;
  // Clear derivative seeds
  val_ind = 0;
  for (i = 0; i < nUnknowns; ++i) {
    // Get variable
    vr = unknowns[i];
    // Get offset, size of variable
    var_off = var_offset[vr];
    var_sz = var_offset[vr + 1] - var_off;
    // Clear elements
    for (j = 0; j < var_sz; ++j) m->d[var_off + j] = 0;
  }
  // Check for evaluation error
  if (flag) return fmi3Error;
  // Successful return
  return fmi3OK;
}
