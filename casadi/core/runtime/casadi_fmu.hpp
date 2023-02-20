// NOLINT(legal/copyright)
typedef struct {
  // Name of instance
  const char* instance_name;
  // Value for all variables
  double v[SZ_MEM];
  // Has any variable been set since last call
  int up_to_date;
  // Buffers for evaluation generated code
  double t;
  double p[N_P];
  double x[N_X];
  double xdot[N_X];
  double u[N_U];
  double y[N_Y];
  // Work vectors for evaluation
  const double* arg[daefun_SZ_ARG];
  double* res[daefun_SZ_RES];
  casadi_int iw[daefun_SZ_IW];
  double w[daefun_SZ_W];
} casadi_fmi_memory;

int evaluate(casadi_fmi_memory* m) {
  // Local variables
  size_t i;
  int mem, flag;
  // Copy states, inputs and parameters to input buffers
  for (i = 0; i < N_X; ++i) m->x[i] = m->v[x_vr[i]];
  for (i = 0; i < N_U; ++i) m->u[i] = m->v[u_vr[i]];
  for (i = 0; i < N_P; ++i) m->p[i] = m->v[p_vr[i]];

  // Evaluate
  m->arg[0] = &m->t;
  m->arg[1] = m->x;
  m->arg[2] = m->u;
  m->arg[3] = 0;
  m->arg[4] = m->p;
  m->res[0] = m->xdot;
  m->res[1] = 0;
  m->res[2] = m->y;
  mem = daefun_checkout();
  flag = daefun(m->arg, m->res, m->iw, m->w, mem);
  daefun_release(mem);

  // Copy from output buffers
  for (i = 0; i < N_X; ++i) m->v[xdot_vr[i]] = m->xdot[i];
  for (i = 0; i < N_Y; ++i) m->v[y_vr[i]] = m->y[i];

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
    daefun_incref();
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
    daefun_decref();
  }
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
