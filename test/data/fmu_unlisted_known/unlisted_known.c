/* Minimal FMI 3.0 Model Exchange FMU with directional and adjoint derivatives:
     der(a) = 2*b,  der(b) = 3*a + 5*b + u
   Like the FMI 3.0 standard, it accepts a known in a derivative only if the
   dependencies list it: time is not listed. */
#include <stdlib.h>
#include "fmi3Functions.h"

enum {VR_TIME, VR_A, VR_B, VR_DA, VR_DB, VR_U, N_VR};

typedef struct {
  double v[N_VR];
} Model;

static void update(Model* m) {
  m->v[VR_DA] = 2*m->v[VR_B];
  m->v[VR_DB] = 3*m->v[VR_A] + 5*m->v[VR_B] + m->v[VR_U];
}

/* Partial derivative of a state derivative with respect to a listed known */
static int partial(fmi3ValueReference unknown, fmi3ValueReference known, double* d) {
  static const double jac[2][3] = {{0, 2, 0}, {3, 5, 1}};
  int i = unknown == VR_DA ? 0 : unknown == VR_DB ? 1 : -1;
  int j = known == VR_A ? 0 : known == VR_B ? 1 : known == VR_U ? 2 : -1;
  if (i < 0 || j < 0) return 1;
  *d = jac[i][j];
  return 0;
}

fmi3Instance fmi3InstantiateModelExchange(fmi3String instanceName,
    fmi3String instantiationToken, fmi3String resourcePath, fmi3Boolean visible,
    fmi3Boolean loggingOn, fmi3InstanceEnvironment instanceEnvironment,
    fmi3LogMessageCallback logMessage) {
  Model* m = (Model*)calloc(1, sizeof(Model));
  m->v[VR_A] = 1;
  m->v[VR_B] = 1;
  return m;
}

void fmi3FreeInstance(fmi3Instance instance) { free(instance); }

fmi3Status fmi3Reset(fmi3Instance instance) { return fmi3Error; }

fmi3Status fmi3EnterInitializationMode(fmi3Instance instance, fmi3Boolean toleranceDefined,
    fmi3Float64 tolerance, fmi3Float64 startTime, fmi3Boolean stopTimeDefined,
    fmi3Float64 stopTime) {
  ((Model*)instance)->v[VR_TIME] = startTime;
  return fmi3OK;
}

fmi3Status fmi3ExitInitializationMode(fmi3Instance instance) { return fmi3OK; }

fmi3Status fmi3EnterContinuousTimeMode(fmi3Instance instance) { return fmi3OK; }

fmi3Status fmi3UpdateDiscreteStates(fmi3Instance instance,
    fmi3Boolean* discreteStatesNeedUpdate, fmi3Boolean* terminateSimulation,
    fmi3Boolean* nominalsOfContinuousStatesChanged, fmi3Boolean* valuesOfContinuousStatesChanged,
    fmi3Boolean* nextEventTimeDefined, fmi3Float64* nextEventTime) {
  *discreteStatesNeedUpdate = fmi3False;
  *terminateSimulation = fmi3False;
  *nominalsOfContinuousStatesChanged = fmi3False;
  *valuesOfContinuousStatesChanged = fmi3False;
  *nextEventTimeDefined = fmi3False;
  return fmi3OK;
}

fmi3Status fmi3SetTime(fmi3Instance instance, fmi3Float64 time) {
  ((Model*)instance)->v[VR_TIME] = time;
  return fmi3OK;
}

fmi3Status fmi3SetFloat64(fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr,
    const fmi3Float64 values[], size_t nValues) {
  Model* m = (Model*)instance;
  size_t k;
  for (k = 0; k < nvr; ++k) {
    if (vr[k] != VR_A && vr[k] != VR_B && vr[k] != VR_U) return fmi3Error;
    m->v[vr[k]] = values[k];
  }
  return fmi3OK;
}

fmi3Status fmi3GetFloat64(fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr,
    fmi3Float64 values[], size_t nValues) {
  Model* m = (Model*)instance;
  size_t k;
  update(m);
  for (k = 0; k < nvr; ++k) {
    if (vr[k] >= N_VR) return fmi3Error;
    values[k] = m->v[vr[k]];
  }
  return fmi3OK;
}

fmi3Status fmi3GetDirectionalDerivative(fmi3Instance instance,
    const fmi3ValueReference unknowns[], size_t nUnknowns,
    const fmi3ValueReference knowns[], size_t nKnowns, const fmi3Float64 seed[], size_t nSeed,
    fmi3Float64 sensitivity[], size_t nSensitivity) {
  size_t i, j;
  double d;
  for (i = 0; i < nUnknowns; ++i) {
    sensitivity[i] = 0;
    for (j = 0; j < nKnowns; ++j) {
      if (partial(unknowns[i], knowns[j], &d)) return fmi3Error;
      sensitivity[i] += d * seed[j];
    }
  }
  return fmi3OK;
}

fmi3Status fmi3GetAdjointDerivative(fmi3Instance instance,
    const fmi3ValueReference unknowns[], size_t nUnknowns,
    const fmi3ValueReference knowns[], size_t nKnowns, const fmi3Float64 seed[], size_t nSeed,
    fmi3Float64 sensitivity[], size_t nSensitivity) {
  size_t i, j;
  double d;
  for (j = 0; j < nKnowns; ++j) {
    sensitivity[j] = 0;
    for (i = 0; i < nUnknowns; ++i) {
      if (partial(unknowns[i], knowns[j], &d)) return fmi3Error;
      sensitivity[j] += d * seed[i];
    }
  }
  return fmi3OK;
}

fmi3Status fmi3SetBoolean(fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr,
    const fmi3Boolean values[], size_t nValues) { return fmi3Error; }

fmi3Status fmi3GetBoolean(fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr,
    fmi3Boolean values[], size_t nValues) { return fmi3Error; }

fmi3Status fmi3SetInt32(fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr,
    const fmi3Int32 values[], size_t nValues) { return fmi3Error; }

fmi3Status fmi3GetInt32(fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr,
    fmi3Int32 values[], size_t nValues) { return fmi3Error; }

fmi3Status fmi3SetString(fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr,
    const fmi3String values[], size_t nValues) { return fmi3Error; }

fmi3Status fmi3GetString(fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr,
    fmi3String values[], size_t nValues) { return fmi3Error; }
