/* Minimal FMI 3.0 Model Exchange FMU with an fmi-ls-dae DAE mode, mimicking an
   OpenModelica --daeMode export:  der(x) = -p*x + z,  0 = z^3 + z - x.
   ODE mode (default): z is solved by Newton internally and der(x) is returned.
   DAE mode (_D_daeMode = true, settable in Configuration Mode only): der(x) and z
   are inputs and the residuals are returned:
     _D_DAEres0 = z^3 + z - x,  _D_DAEres1 = 2*(der(x) + p*x - z) */
#include <stdlib.h>
#include "fmi3Functions.h"

enum {VR_TIME, VR_X, VR_DX, VR_Z, VR_P, VR_DAEMODE, VR_RES0, VR_RES1, N_VR};
enum {MODE_INSTANTIATED, MODE_CONFIGURATION, MODE_INITIALIZATION, MODE_CONTINUOUS};

typedef struct {
  double v[N_VR];
  int dae_mode, mode;
} Model;

static void update(Model* m) {
  double x = m->v[VR_X], p = m->v[VR_P], z = m->v[VR_Z];
  if (!m->dae_mode) {
    int k;
    for (k = 0; k < 50; ++k) z -= (z*z*z + z - x) / (3*z*z + 1);
    m->v[VR_Z] = z;
    m->v[VR_DX] = -p*x + z;
  }
  m->v[VR_RES0] = z*z*z + z - x;
  m->v[VR_RES1] = 2*(m->v[VR_DX] + p*x - z);
}

fmi3Instance fmi3InstantiateModelExchange(fmi3String instanceName,
    fmi3String instantiationToken, fmi3String resourcePath, fmi3Boolean visible,
    fmi3Boolean loggingOn, fmi3InstanceEnvironment instanceEnvironment,
    fmi3LogMessageCallback logMessage) {
  Model* m = (Model*)calloc(1, sizeof(Model));
  m->v[VR_X] = 1;
  m->v[VR_P] = 2;
  m->mode = MODE_INSTANTIATED;
  return m;
}

void fmi3FreeInstance(fmi3Instance instance) { free(instance); }

fmi3Status fmi3Reset(fmi3Instance instance) { return fmi3Error; }

fmi3Status fmi3EnterConfigurationMode(fmi3Instance instance) {
  Model* m = (Model*)instance;
  if (m->mode != MODE_INSTANTIATED) return fmi3Error;
  m->mode = MODE_CONFIGURATION;
  return fmi3OK;
}

fmi3Status fmi3ExitConfigurationMode(fmi3Instance instance) {
  Model* m = (Model*)instance;
  if (m->mode != MODE_CONFIGURATION) return fmi3Error;
  m->mode = MODE_INSTANTIATED;
  return fmi3OK;
}

fmi3Status fmi3EnterInitializationMode(fmi3Instance instance, fmi3Boolean toleranceDefined,
    fmi3Float64 tolerance, fmi3Float64 startTime, fmi3Boolean stopTimeDefined,
    fmi3Float64 stopTime) {
  Model* m = (Model*)instance;
  m->v[VR_TIME] = startTime;
  m->mode = MODE_INITIALIZATION;
  return fmi3OK;
}

fmi3Status fmi3ExitInitializationMode(fmi3Instance instance) { return fmi3OK; }

fmi3Status fmi3EnterContinuousTimeMode(fmi3Instance instance) {
  ((Model*)instance)->mode = MODE_CONTINUOUS;
  return fmi3OK;
}

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
    switch (vr[k]) {
      case VR_X: case VR_P: break;
      case VR_DX: case VR_Z: if (m->dae_mode) break; /* fall-through */
      default: return fmi3Error;
    }
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
    if (vr[k] >= N_VR || vr[k] == VR_DAEMODE) return fmi3Error;
    values[k] = m->v[vr[k]];
  }
  return fmi3OK;
}

fmi3Status fmi3SetBoolean(fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr,
    const fmi3Boolean values[], size_t nValues) {
  Model* m = (Model*)instance;
  if (nvr != 1 || vr[0] != VR_DAEMODE || m->mode != MODE_CONFIGURATION) return fmi3Error;
  m->dae_mode = values[0];
  return fmi3OK;
}

fmi3Status fmi3GetBoolean(fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr,
    fmi3Boolean values[], size_t nValues) {
  Model* m = (Model*)instance;
  if (nvr != 1 || vr[0] != VR_DAEMODE) return fmi3Error;
  values[0] = m->dae_mode;
  return fmi3OK;
}

fmi3Status fmi3SetInt32(fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr,
    const fmi3Int32 values[], size_t nValues) { return fmi3Error; }

fmi3Status fmi3GetInt32(fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr,
    fmi3Int32 values[], size_t nValues) { return fmi3Error; }

fmi3Status fmi3SetString(fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr,
    const fmi3String values[], size_t nValues) { return fmi3Error; }

fmi3Status fmi3GetString(fmi3Instance instance, const fmi3ValueReference vr[], size_t nvr,
    fmi3String values[], size_t nValues) { return fmi3Error; }
