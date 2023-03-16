#ifndef fmi3FunctionTypes_h
#define fmi3FunctionTypes_h

#include "fmi3PlatformTypes.h"

/*
This header file defines the data and function types of FMI 3.0.
It must be used when compiling an FMU or an FMI importer.

Copyright (C) 2011 MODELISAR consortium,
              2012-2022 Modelica Association Project "FMI"
              All rights reserved.

This file is licensed by the copyright holders under the 2-Clause BSD License
(https://opensource.org/licenses/BSD-2-Clause):

----------------------------------------------------------------------------
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
----------------------------------------------------------------------------
*/

#ifdef __cplusplus
extern "C" {
#endif

/* Include stddef.h, in order that size_t etc. is defined */
#include <stddef.h>


/* Type definitions */

/* tag::Status[] */
typedef enum {
    fmi3OK,
    fmi3Warning,
    fmi3Discard,
    fmi3Error,
    fmi3Fatal,
} fmi3Status;
/* end::Status[] */

/* tag::DependencyKind[] */
typedef enum {
    fmi3Independent,
    fmi3Constant,
    fmi3Fixed,
    fmi3Tunable,
    fmi3Discrete,
    fmi3Dependent
} fmi3DependencyKind;
/* end::DependencyKind[] */

/* tag::IntervalQualifier[] */
typedef enum {
    fmi3IntervalNotYetKnown,
    fmi3IntervalUnchanged,
    fmi3IntervalChanged
} fmi3IntervalQualifier;
/* end::IntervalQualifier[] */

/* tag::CallbackLogMessage[] */
typedef void  (*fmi3LogMessageCallback) (fmi3InstanceEnvironment instanceEnvironment,
                                         fmi3Status status,
                                         fmi3String category,
                                         fmi3String message);
/* end::CallbackLogMessage[] */

/* tag::CallbackClockUpdate[] */
typedef void (*fmi3ClockUpdateCallback) (
    fmi3InstanceEnvironment  instanceEnvironment);
/* end::CallbackClockUpdate[] */

/* tag::CallbackIntermediateUpdate[] */
typedef void (*fmi3IntermediateUpdateCallback) (
    fmi3InstanceEnvironment instanceEnvironment,
    fmi3Float64  intermediateUpdateTime,
    fmi3Boolean  intermediateVariableSetRequested,
    fmi3Boolean  intermediateVariableGetAllowed,
    fmi3Boolean  intermediateStepFinished,
    fmi3Boolean  canReturnEarly,
    fmi3Boolean* earlyReturnRequested,
    fmi3Float64* earlyReturnTime);
/* end::CallbackIntermediateUpdate[] */

/* tag::CallbackPreemptionLock[] */
typedef void (*fmi3LockPreemptionCallback)   ();
typedef void (*fmi3UnlockPreemptionCallback) ();
/* end::CallbackPreemptionLock[] */

/* Define fmi3 function pointer types to simplify dynamic loading */

/***************************************************
Types for Common Functions
****************************************************/

/* Inquire version numbers and setting logging status */
/* tag::GetVersion[] */
typedef const char* fmi3GetVersionTYPE(void);
/* end::GetVersion[] */

/* tag::SetDebugLogging[] */
typedef fmi3Status fmi3SetDebugLoggingTYPE(fmi3Instance instance,
                                           fmi3Boolean loggingOn,
                                           size_t nCategories,
                                           const fmi3String categories[]);
/* end::SetDebugLogging[] */

/* Creation and destruction of FMU instances and setting debug status */
/* tag::Instantiate[] */
typedef fmi3Instance fmi3InstantiateModelExchangeTYPE(
    fmi3String                 instanceName,
    fmi3String                 instantiationToken,
    fmi3String                 resourcePath,
    fmi3Boolean                visible,
    fmi3Boolean                loggingOn,
    fmi3InstanceEnvironment    instanceEnvironment,
    fmi3LogMessageCallback     logMessage);

typedef fmi3Instance fmi3InstantiateCoSimulationTYPE(
    fmi3String                     instanceName,
    fmi3String                     instantiationToken,
    fmi3String                     resourcePath,
    fmi3Boolean                    visible,
    fmi3Boolean                    loggingOn,
    fmi3Boolean                    eventModeUsed,
    fmi3Boolean                    earlyReturnAllowed,
    const fmi3ValueReference       requiredIntermediateVariables[],
    size_t                         nRequiredIntermediateVariables,
    fmi3InstanceEnvironment        instanceEnvironment,
    fmi3LogMessageCallback         logMessage,
    fmi3IntermediateUpdateCallback intermediateUpdate);

typedef fmi3Instance fmi3InstantiateScheduledExecutionTYPE(
    fmi3String                     instanceName,
    fmi3String                     instantiationToken,
    fmi3String                     resourcePath,
    fmi3Boolean                    visible,
    fmi3Boolean                    loggingOn,
    fmi3InstanceEnvironment        instanceEnvironment,
    fmi3LogMessageCallback         logMessage,
    fmi3ClockUpdateCallback        clockUpdate,
    fmi3LockPreemptionCallback     lockPreemption,
    fmi3UnlockPreemptionCallback   unlockPreemption);
/* end::Instantiate[] */

/* tag::FreeInstance[] */
typedef void fmi3FreeInstanceTYPE(fmi3Instance instance);
/* end::FreeInstance[] */

/* Enter and exit initialization mode, enter event mode, terminate and reset */
/* tag::EnterInitializationMode[] */
typedef fmi3Status fmi3EnterInitializationModeTYPE(fmi3Instance instance,
                                                   fmi3Boolean toleranceDefined,
                                                   fmi3Float64 tolerance,
                                                   fmi3Float64 startTime,
                                                   fmi3Boolean stopTimeDefined,
                                                   fmi3Float64 stopTime);
/* end::EnterInitializationMode[] */

/* tag::ExitInitializationMode[] */
typedef fmi3Status fmi3ExitInitializationModeTYPE(fmi3Instance instance);
/* end::ExitInitializationMode[] */

/* tag::EnterEventMode[] */
typedef fmi3Status fmi3EnterEventModeTYPE(fmi3Instance instance);
/* end::EnterEventMode[] */

/* tag::Terminate[] */
typedef fmi3Status fmi3TerminateTYPE(fmi3Instance instance);
/* end::Terminate[] */

/* tag::Reset[] */
typedef fmi3Status fmi3ResetTYPE(fmi3Instance instance);
/* end::Reset[] */

/* Getting and setting variable values */
/* tag::Getters[] */
typedef fmi3Status fmi3GetFloat32TYPE(fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Float32 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetFloat64TYPE(fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Float64 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetInt8TYPE   (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Int8 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetUInt8TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3UInt8 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetInt16TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Int16 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetUInt16TYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3UInt16 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetInt32TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Int32 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetUInt32TYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3UInt32 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetInt64TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Int64 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetUInt64TYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3UInt64 values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetBooleanTYPE(fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Boolean values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetStringTYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3String values[],
                                      size_t nValues);

typedef fmi3Status fmi3GetBinaryTYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      size_t valueSizes[],
                                      fmi3Binary values[],
                                      size_t nValues);
/* end::Getters[] */

/* tag::GetClock[] */
typedef fmi3Status fmi3GetClockTYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      fmi3Clock values[]);
/* end::GetClock[] */

/* tag::Setters[] */
typedef fmi3Status fmi3SetFloat32TYPE(fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Float32 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetFloat64TYPE(fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Float64 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetInt8TYPE   (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Int8 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetUInt8TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3UInt8 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetInt16TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Int16 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetUInt16TYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3UInt16 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetInt32TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Int32 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetUInt32TYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3UInt32 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetInt64TYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Int64 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetUInt64TYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3UInt64 values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetBooleanTYPE(fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Boolean values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetStringTYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3String values[],
                                      size_t nValues);

typedef fmi3Status fmi3SetBinaryTYPE (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const size_t valueSizes[],
                                      const fmi3Binary values[],
                                      size_t nValues);
/* end::Setters[] */
/* tag::SetClock[] */
typedef fmi3Status fmi3SetClockTYPE  (fmi3Instance instance,
                                      const fmi3ValueReference valueReferences[],
                                      size_t nValueReferences,
                                      const fmi3Clock values[]);
/* end::SetClock[] */

/* Getting Variable Dependency Information */
/* tag::GetNumberOfVariableDependencies[] */
typedef fmi3Status fmi3GetNumberOfVariableDependenciesTYPE(fmi3Instance instance,
                                                           fmi3ValueReference valueReference,
                                                           size_t* nDependencies);
/* end::GetNumberOfVariableDependencies[] */

/* tag::GetVariableDependencies[] */
typedef fmi3Status fmi3GetVariableDependenciesTYPE(fmi3Instance instance,
                                                   fmi3ValueReference dependent,
                                                   size_t elementIndicesOfDependent[],
                                                   fmi3ValueReference independents[],
                                                   size_t elementIndicesOfIndependents[],
                                                   fmi3DependencyKind dependencyKinds[],
                                                   size_t nDependencies);
/* end::GetVariableDependencies[] */

/* Getting and setting the internal FMU state */
/* tag::GetFMUState[] */
typedef fmi3Status fmi3GetFMUStateTYPE (fmi3Instance instance, fmi3FMUState* FMUState);
/* end::GetFMUState[] */

/* tag::SetFMUState[] */
typedef fmi3Status fmi3SetFMUStateTYPE (fmi3Instance instance, fmi3FMUState  FMUState);
/* end::SetFMUState[] */

/* tag::FreeFMUState[] */
typedef fmi3Status fmi3FreeFMUStateTYPE(fmi3Instance instance, fmi3FMUState* FMUState);
/* end::FreeFMUState[] */

/* tag::SerializedFMUStateSize[] */
typedef fmi3Status fmi3SerializedFMUStateSizeTYPE(fmi3Instance instance,
                                                  fmi3FMUState FMUState,
                                                  size_t* size);
/* end::SerializedFMUStateSize[] */

/* tag::SerializeFMUState[] */
typedef fmi3Status fmi3SerializeFMUStateTYPE     (fmi3Instance instance,
                                                  fmi3FMUState FMUState,
                                                  fmi3Byte serializedState[],
                                                  size_t size);
/* end::SerializeFMUState[] */

/* tag::DeserializeFMUState[] */
typedef fmi3Status fmi3DeserializeFMUStateTYPE   (fmi3Instance instance,
                                                  const fmi3Byte serializedState[],
                                                  size_t size,
                                                  fmi3FMUState* FMUState);
/* end::DeserializeFMUState[] */

/* Getting partial derivatives */
/* tag::GetDirectionalDerivative[] */
typedef fmi3Status fmi3GetDirectionalDerivativeTYPE(fmi3Instance instance,
                                                    const fmi3ValueReference unknowns[],
                                                    size_t nUnknowns,
                                                    const fmi3ValueReference knowns[],
                                                    size_t nKnowns,
                                                    const fmi3Float64 seed[],
                                                    size_t nSeed,
                                                    fmi3Float64 sensitivity[],
                                                    size_t nSensitivity);
/* end::GetDirectionalDerivative[] */

/* tag::GetAdjointDerivative[] */
typedef fmi3Status fmi3GetAdjointDerivativeTYPE(fmi3Instance instance,
                                                const fmi3ValueReference unknowns[],
                                                size_t nUnknowns,
                                                const fmi3ValueReference knowns[],
                                                size_t nKnowns,
                                                const fmi3Float64 seed[],
                                                size_t nSeed,
                                                fmi3Float64 sensitivity[],
                                                size_t nSensitivity);
/* end::GetAdjointDerivative[] */

/* Entering and exiting the Configuration or Reconfiguration Mode */

/* tag::EnterConfigurationMode[] */
typedef fmi3Status fmi3EnterConfigurationModeTYPE(fmi3Instance instance);
/* end::EnterConfigurationMode[] */

/* tag::ExitConfigurationMode[] */
typedef fmi3Status fmi3ExitConfigurationModeTYPE(fmi3Instance instance);
/* end::ExitConfigurationMode[] */

/* tag::GetIntervalDecimal[] */
typedef fmi3Status fmi3GetIntervalDecimalTYPE(fmi3Instance instance,
                                              const fmi3ValueReference valueReferences[],
                                              size_t nValueReferences,
                                              fmi3Float64 intervals[],
                                              fmi3IntervalQualifier qualifiers[]);
/* end::GetIntervalDecimal[] */

/* tag::GetIntervalFraction[] */
typedef fmi3Status fmi3GetIntervalFractionTYPE(fmi3Instance instance,
                                               const fmi3ValueReference valueReferences[],
                                               size_t nValueReferences,
                                               fmi3UInt64 counters[],
                                               fmi3UInt64 resolutions[],
                                               fmi3IntervalQualifier qualifiers[]);
/* end::GetIntervalFraction[] */

/* tag::GetShiftDecimal[] */
typedef fmi3Status fmi3GetShiftDecimalTYPE(fmi3Instance instance,
                                           const fmi3ValueReference valueReferences[],
                                           size_t nValueReferences,
                                           fmi3Float64 shifts[]);
/* end::GetShiftDecimal[] */

/* tag::GetShiftFraction[] */
typedef fmi3Status fmi3GetShiftFractionTYPE(fmi3Instance instance,
                                            const fmi3ValueReference valueReferences[],
                                            size_t nValueReferences,
                                            fmi3UInt64 counters[],
                                            fmi3UInt64 resolutions[]);
/* end::GetShiftFraction[] */

/* tag::SetIntervalDecimal[] */
typedef fmi3Status fmi3SetIntervalDecimalTYPE(fmi3Instance instance,
                                              const fmi3ValueReference valueReferences[],
                                              size_t nValueReferences,
                                              const fmi3Float64 intervals[]);
/* end::SetIntervalDecimal[] */

/* tag::SetIntervalFraction[] */
typedef fmi3Status fmi3SetIntervalFractionTYPE(fmi3Instance instance,
                                               const fmi3ValueReference valueReferences[],
                                               size_t nValueReferences,
                                               const fmi3UInt64 counters[],
                                               const fmi3UInt64 resolutions[]);
/* end::SetIntervalFraction[] */

/* tag::SetShiftDecimal[] */
typedef fmi3Status fmi3SetShiftDecimalTYPE(fmi3Instance instance,
                                           const fmi3ValueReference valueReferences[],
                                           size_t nValueReferences,
                                           const fmi3Float64 shifts[]);
/* end::SetShiftDecimal[] */

/* tag::SetShiftFraction[] */
typedef fmi3Status fmi3SetShiftFractionTYPE(fmi3Instance instance,
                                            const fmi3ValueReference valueReferences[],
                                            size_t nValueReferences,
                                            const fmi3UInt64 counters[],
                                            const fmi3UInt64 resolutions[]);
/* end::SetShiftFraction[] */

/* tag::EvaluateDiscreteStates[] */
typedef fmi3Status fmi3EvaluateDiscreteStatesTYPE(fmi3Instance instance);
/* end::EvaluateDiscreteStates[] */

/* tag::UpdateDiscreteStates[] */
typedef fmi3Status fmi3UpdateDiscreteStatesTYPE(fmi3Instance instance,
                                                fmi3Boolean* discreteStatesNeedUpdate,
                                                fmi3Boolean* terminateSimulation,
                                                fmi3Boolean* nominalsOfContinuousStatesChanged,
                                                fmi3Boolean* valuesOfContinuousStatesChanged,
                                                fmi3Boolean* nextEventTimeDefined,
                                                fmi3Float64* nextEventTime);
/* end::UpdateDiscreteStates[] */

/***************************************************
Types for Functions for Model Exchange
****************************************************/

/* tag::EnterContinuousTimeMode[] */
typedef fmi3Status fmi3EnterContinuousTimeModeTYPE(fmi3Instance instance);
/* end::EnterContinuousTimeMode[] */

/* tag::CompletedIntegratorStep[] */
typedef fmi3Status fmi3CompletedIntegratorStepTYPE(fmi3Instance instance,
                                                   fmi3Boolean  noSetFMUStatePriorToCurrentPoint,
                                                   fmi3Boolean* enterEventMode,
                                                   fmi3Boolean* terminateSimulation);
/* end::CompletedIntegratorStep[] */

/* Providing independent variables and re-initialization of caching */
/* tag::SetTime[] */
typedef fmi3Status fmi3SetTimeTYPE(fmi3Instance instance, fmi3Float64 time);
/* end::SetTime[] */

/* tag::SetContinuousStates[] */
typedef fmi3Status fmi3SetContinuousStatesTYPE(fmi3Instance instance,
                                               const fmi3Float64 continuousStates[],
                                               size_t nContinuousStates);
/* end::SetContinuousStates[] */

/* Evaluation of the model equations */
/* tag::GetDerivatives[] */
typedef fmi3Status fmi3GetContinuousStateDerivativesTYPE(fmi3Instance instance,
                                                         fmi3Float64 derivatives[],
                                                         size_t nContinuousStates);
/* end::GetDerivatives[] */

/* tag::GetEventIndicators[] */
typedef fmi3Status fmi3GetEventIndicatorsTYPE(fmi3Instance instance,
                                              fmi3Float64 eventIndicators[],
                                              size_t nEventIndicators);
/* end::GetEventIndicators[] */

/* tag::GetContinuousStates[] */
typedef fmi3Status fmi3GetContinuousStatesTYPE(fmi3Instance instance,
                                               fmi3Float64 continuousStates[],
                                               size_t nContinuousStates);
/* end::GetContinuousStates[] */

/* tag::GetNominalsOfContinuousStates[] */
typedef fmi3Status fmi3GetNominalsOfContinuousStatesTYPE(fmi3Instance instance,
                                                         fmi3Float64 nominals[],
                                                         size_t nContinuousStates);
/* end::GetNominalsOfContinuousStates[] */

/* tag::GetNumberOfEventIndicators[] */
typedef fmi3Status fmi3GetNumberOfEventIndicatorsTYPE(fmi3Instance instance,
                                                      size_t* nEventIndicators);
/* end::GetNumberOfEventIndicators[] */

/* tag::GetNumberOfContinuousStates[] */
typedef fmi3Status fmi3GetNumberOfContinuousStatesTYPE(fmi3Instance instance,
                                                       size_t* nContinuousStates);
/* end::GetNumberOfContinuousStates[] */

/***************************************************
Types for Functions for Co-Simulation
****************************************************/

/* Simulating the FMU */

/* tag::EnterStepMode[] */
typedef fmi3Status fmi3EnterStepModeTYPE(fmi3Instance instance);
/* end::EnterStepMode[] */

/* tag::GetOutputDerivatives[] */
typedef fmi3Status fmi3GetOutputDerivativesTYPE(fmi3Instance instance,
                                                const fmi3ValueReference valueReferences[],
                                                size_t nValueReferences,
                                                const fmi3Int32 orders[],
                                                fmi3Float64 values[],
                                                size_t nValues);
/* end::GetOutputDerivatives[] */

/* tag::DoStep[] */
typedef fmi3Status fmi3DoStepTYPE(fmi3Instance instance,
                                  fmi3Float64 currentCommunicationPoint,
                                  fmi3Float64 communicationStepSize,
                                  fmi3Boolean noSetFMUStatePriorToCurrentPoint,
                                  fmi3Boolean* eventHandlingNeeded,
                                  fmi3Boolean* terminateSimulation,
                                  fmi3Boolean* earlyReturn,
                                  fmi3Float64* lastSuccessfulTime);
/* end::DoStep[] */

/***************************************************
Types for Functions for Scheduled Execution
****************************************************/

/* tag::ActivateModelPartition[] */
typedef fmi3Status fmi3ActivateModelPartitionTYPE(fmi3Instance instance,
                                                  fmi3ValueReference clockReference,
                                                  fmi3Float64 activationTime);
/* end::ActivateModelPartition[] */

#ifdef __cplusplus
}  /* end of extern "C" { */
#endif

#endif /* fmi3FunctionTypes_h */
