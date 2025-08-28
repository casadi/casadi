# Eigen library for vectors and matrices
find_package(Eigen3 REQUIRED)
# Guanaqo library for all kinds of utilities
find_package(guanaqo REQUIRED)

function(alpaqa_configure_visibility target)
    set_target_properties(${target} PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                                               C_VISIBILITY_PRESET "hidden"
                                               VISIBILITY_INLINES_HIDDEN true)
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
        target_link_options(${target} PRIVATE
            $<$<LINK_LANGUAGE:C,CXX>:LINKER:--exclude-libs,ALL>)
    endif()
endfunction()

add_library(alpaqa
    "alpaqa/src/accelerators/lbfgs.cpp"
    "alpaqa/src/problem/problem-counters.cpp"
    "alpaqa/src/problem/ocproblem-counters.cpp"
    "alpaqa/src/problem/type-erased-problem.cpp"
    "alpaqa/src/outer/alm.cpp"
    "alpaqa/src/outer/internal/alm-helpers.cpp"
    "alpaqa/src/inner/panoc.cpp"
    "alpaqa/src/inner/fista.cpp"
    "alpaqa/src/inner/zerofpr.cpp"
    "alpaqa/src/inner/pantr.cpp"
    "alpaqa/src/panoc-alm.cpp"
    "alpaqa/src/panoc-convex-newton-alm.cpp"
    "alpaqa/src/panoc-anderson-alm.cpp"
    "alpaqa/src/structured-panoc-alm.cpp"
    "alpaqa/src/structured-zerofpr-alm.cpp"
    "alpaqa/src/zerofpr-alm.cpp"
    "alpaqa/src/zerofpr-convex-newton-alm.cpp"
    "alpaqa/src/zerofpr-anderson-alm.cpp"
    "alpaqa/src/newton-tr-pantr-alm.cpp"
    "alpaqa/src/inner/internal/solverstatus.cpp"
    "alpaqa/src/inner/directions/panoc/structured-lbfgs.cpp"
    "alpaqa/src/inner/directions/panoc/structured-newton.cpp"
    "alpaqa/src/inner/directions/panoc/convex-newton.cpp"
    "alpaqa/src/inner/internal/panoc-helpers.cpp"
    "alpaqa/src/params/params.cpp"
)
if (ALPAQA_WITH_OCP)
    target_sources(alpaqa PRIVATE
        "alpaqa/src/inner/panoc-ocp.cpp"
        "alpaqa/src/problem/ocproblem.cpp")
endif()
if (ALPAQA_WITH_JSON)
    find_package(nlohmann_json REQUIRED)
    target_sources(alpaqa PRIVATE
        "alpaqa/src/params/json.cpp")
    target_link_libraries(alpaqa
        PUBLIC $<BUILD_INTERFACE:nlohmann_json::nlohmann_json>
        # Disable in the install interface, enabled again in 
        # cmake/CoreConfig.cmake.in
        PRIVATE $<INSTALL_INTERFACE:nlohmann_json::nlohmann_json>)
endif()

target_compile_features(alpaqa PUBLIC cxx_std_20)
if (ALPAQA_WITH_CXX_23)
    # Use C++23 to build everything, but don't require C++23 for installed
    # versions of the library (see cmake/CoreConfig.cmake.in)
    target_compile_features(alpaqa PUBLIC $<BUILD_INTERFACE:cxx_std_23>)
endif()
target_include_directories(alpaqa PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/alpaqa/include>
    $<INSTALL_INTERFACE:${ALPAQA_INSTALL_INCLUDEDIR}>)
if (ALPAQA_DEBUG_CHECKS_EIGEN)
    target_compile_definitions(alpaqa PUBLIC
        $<$<CONFIG:Debug>:EIGEN_INITIALIZE_MATRICES_BY_NAN>
        $<$<CONFIG:Debug>:EIGEN_RUNTIME_NO_MALLOC>)
endif()
target_link_libraries(alpaqa PUBLIC Eigen3::Eigen guanaqo::guanaqo guanaqo::linalg)
target_link_libraries(alpaqa PRIVATE warnings)
alpaqa_configure_visibility(alpaqa)
target_compile_definitions(alpaqa PUBLIC
    $<$<BOOL:${ALPAQA_WITH_QUAD_PRECISION}>:ALPAQA_WITH_QUAD_PRECISION>
    $<$<BOOL:${ALPAQA_WITH_SINGLE_PRECISION}>:ALPAQA_WITH_SINGLE_PRECISION>
    $<$<BOOL:${ALPAQA_WITH_LONG_DOUBLE}>:ALPAQA_WITH_LONG_DOUBLE>
    $<$<BOOL:${ALPAQA_NO_DLCLOSE}>:ALPAQA_NO_DLCLOSE>
    $<$<BOOL:${ALPAQA_WITH_JSON}>:ALPAQA_WITH_JSON>
    $<$<BOOL:${ALPAQA_WITH_OCP}>:ALPAQA_WITH_OCP>
)
target_link_libraries(alpaqa PUBLIC
    $<$<BOOL:${ALPAQA_WITH_QUAD_PRECISION}>:quadmath>)
if (MSVC)
    target_compile_options(alpaqa PUBLIC "/utf-8")
endif()
if (ALPAQA_WITH_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
    target_link_libraries(alpaqa PUBLIC OpenMP::OpenMP_CXX)
endif()
if (ALPAQA_DONT_PARALLELIZE_EIGEN)
    target_compile_definitions(alpaqa PUBLIC EIGEN_DONT_PARALLELIZE)
endif()
if (ALPAQA_WITH_BLAS)
    find_package(BLAS REQUIRED)
    target_link_libraries(alpaqa PUBLIC BLAS::BLAS)
    target_compile_definitions(alpaqa PUBLIC EIGEN_USE_BLAS)
endif()
set_target_properties(alpaqa PROPERTIES SOVERSION ${PROJECT_VERSION})
add_library(alpaqa::alpaqa ALIAS alpaqa)
list(APPEND ALPAQA_INSTALL_TARGETS alpaqa)

# DLL import/export
include(GenerateExportHeader)
generate_export_header(alpaqa
    EXPORT_FILE_NAME export/alpaqa/export.h)
target_include_directories(alpaqa PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)

# Build time
if (ALPAQA_WITH_ACCURATE_BUILD_TIME)
    set(ALPAQA_BUILD_TIME_CPP "${CMAKE_CURRENT_BINARY_DIR}/alpaqa-build-time.cpp")
    add_custom_target(alpaqa-build-time-generate
        BYPRODUCTS ${ALPAQA_BUILD_TIME_CPP}
        COMMAND ${CMAKE_COMMAND}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BuildTime.cmake
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    target_sources(alpaqa PRIVATE ${ALPAQA_BUILD_TIME_CPP})
    add_dependencies(alpaqa alpaqa-build-time-generate)
else()
    set(ALPAQA_BUILD_TIME_CPP "${CMAKE_CURRENT_BINARY_DIR}/alpaqa-build-time.cpp")
    if (NOT EXISTS ${ALPAQA_BUILD_TIME_CPP})
        include(cmake/BuildTime.cmake)
    endif()
    target_sources(alpaqa PRIVATE ${ALPAQA_BUILD_TIME_CPP})
endif()

# CasADi
if (ALPAQA_WITH_EXTERNAL_CASADI)
    find_package(casadi REQUIRED)
    find_package(Threads REQUIRED)
    if (TARGET casadi AND NOT TARGET casadi::casadi)
        add_library(casadi::casadi ALIAS casadi)  # Support CasADi <3.7
    endif()
endif()
if (ALPAQA_WITH_CASADI)
    # Normal NLPs
    add_library(casadi-loader
        "interop/casadi/src/CasADiProblem.cpp"
        "interop/casadi/include/alpaqa/casadi/CasADiFunctionWrapper.hpp")
    target_include_directories(casadi-loader PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/casadi/include>
        $<INSTALL_INTERFACE:${ALPAQA_INSTALL_INCLUDEDIR}>)
    target_compile_definitions(casadi-loader PUBLIC
        $<$<BOOL:${ALPAQA_WITH_CASADI}>:ALPAQA_WITH_CASADI>
        $<$<BOOL:${ALPAQA_WITH_EXTERNAL_CASADI}>:ALPAQA_WITH_EXTERNAL_CASADI>)
    set_property(TARGET casadi-loader
        PROPERTY OUTPUT_NAME alpaqa-casadi-loader)
    target_link_libraries(casadi-loader
        PUBLIC  alpaqa::alpaqa
        PRIVATE alpaqa::warnings)
    if (ALPAQA_WITH_EXTERNAL_CASADI)
        target_sources(casadi-loader PRIVATE
            "interop/casadi/src/CompleteCasADiProblem.cpp")
        target_link_libraries(casadi-loader
            PRIVATE casadi::casadi Threads::Threads)
    else()
        target_sources(casadi-loader PRIVATE
            "interop/casadi/src/casadi-external-function.cpp"
            "interop/casadi/src/casadi-external-function.cpp")
    endif()
    alpaqa_configure_visibility(casadi-loader)
    set_target_properties(casadi-loader PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::casadi-loader ALIAS casadi-loader)
    list(APPEND ALPAQA_INSTALL_TARGETS casadi-loader)
    generate_export_header(casadi-loader
        EXPORT_FILE_NAME export/alpaqa/casadi-loader-export.h)
    target_include_directories(casadi-loader PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)

    message(STATUS "Compiling WITH CasADi support")
else()
    message(STATUS "Compiling WITHOUT CasADi support")
endif()

# Optimal control problems
if (ALPAQA_WITH_OCP AND ALPAQA_WITH_CASADI_OCP)
    add_library(casadi-ocp-loader "interop/casadi/src/CasADiControlProblem.cpp")
    target_compile_definitions(casadi-loader PUBLIC
        $<$<BOOL:${ALPAQA_WITH_CASADI_OCP}>:ALPAQA_WITH_CASADI_OCP>)
    set_property(TARGET casadi-ocp-loader
        PROPERTY OUTPUT_NAME alpaqa-casadi-ocp-loader)
    target_link_libraries(casadi-ocp-loader
        PUBLIC  alpaqa::alpaqa alpaqa::casadi-loader
        PRIVATE alpaqa::warnings)
    if (ALPAQA_WITH_EXTERNAL_CASADI)
        target_link_libraries(casadi-ocp-loader
            PRIVATE casadi::casadi Threads::Threads)
    endif()
    alpaqa_configure_visibility(casadi-ocp-loader)
    set_target_properties(casadi-ocp-loader PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::casadi-ocp-loader ALIAS casadi-ocp-loader)
    list(APPEND ALPAQA_INSTALL_TARGETS casadi-ocp-loader)
    generate_export_header(casadi-ocp-loader
        EXPORT_FILE_NAME export/alpaqa/casadi-ocp-loader-export.h)
    target_include_directories(casadi-ocp-loader PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)
endif()

add_library(dl-api INTERFACE)
target_include_directories(dl-api INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/dl-api/include>
    $<INSTALL_INTERFACE:${ALPAQA_INSTALL_INCLUDEDIR}>)
if (MSVC)
    target_compile_options(dl-api INTERFACE "/utf-8")
endif()
add_library(alpaqa::dl-api ALIAS dl-api)
include(cmake/dl-problem.cmake)
list(APPEND ALPAQA_INSTALL_TARGETS dl-api)

if (ALPAQA_WITH_DL)
    add_library(dl-loader "interop/dl/src/dl-problem.cpp")
    target_include_directories(dl-loader PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/dl/include>
        $<INSTALL_INTERFACE:${ALPAQA_INSTALL_INCLUDEDIR}>)
    target_compile_definitions(dl-loader PUBLIC ALPAQA_WITH_DL)
    set_property(TARGET dl-loader
        PROPERTY OUTPUT_NAME alpaqa-dl-loader)
    target_link_libraries(dl-loader
        PUBLIC alpaqa::alpaqa alpaqa::dl-api
        PRIVATE alpaqa::warnings)
    alpaqa_configure_visibility(dl-loader)
    set_target_properties(dl-loader PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::dl-loader ALIAS dl-loader)
    list(APPEND ALPAQA_INSTALL_TARGETS dl-loader)
    generate_export_header(dl-loader
        EXPORT_FILE_NAME export/alpaqa/dl-loader-export.h)
    target_include_directories(dl-loader PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)
endif()

if (ALPAQA_WITH_CUTEST AND CMAKE_DL_LIBS)
    add_library(cutest-interface "interop/cutest/src/cutest-loader.cpp")
    target_include_directories(cutest-interface PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/cutest/include>
        $<INSTALL_INTERFACE:${ALPAQA_INSTALL_INCLUDEDIR}>)
    set_property(TARGET cutest-interface
        PROPERTY OUTPUT_NAME alpaqa-cutest-interface)
    target_link_libraries(cutest-interface
        PUBLIC alpaqa::alpaqa
        PRIVATE ${CMAKE_DL_LIBS} alpaqa::warnings)
    if (NOT ALPAQA_WITH_CXX_23)
        message(FATAL_ERROR "ALPAQA_WITH_CUTEST requires C++23 support")
    endif()
    target_compile_definitions(cutest-interface PUBLIC ALPAQA_WITH_CUTEST)
    alpaqa_configure_visibility(cutest-interface)
    set_target_properties(cutest-interface PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::cutest-interface ALIAS cutest-interface)
    list(APPEND ALPAQA_INSTALL_TARGETS cutest-interface)
    generate_export_header(cutest-interface
        EXPORT_FILE_NAME export/alpaqa/cutest-interface-export.h)
    target_include_directories(cutest-interface PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)
    message(STATUS "Compiling WITH CUTEst support")
elseif(ALPAQA_WITH_CUTEST)
    message(STATUS "Compiling WITHOUT CUTEst support: missing dlopen")
else()
    message(STATUS "Compiling WITHOUT CUTEst support")
endif()

if (ALPAQA_WITH_IPOPT)
    find_package(Ipopt REQUIRED)
    message(STATUS "Compiling Ipopt problem adapter")
    add_library(ipopt-adapter "interop/ipopt/src/ipopt-adapter.cpp")
    if (ALPAQA_WITH_DRIVERS)
        target_sources(ipopt-adapter PRIVATE
            "interop/ipopt/src/ipopt-params.cpp")
    endif()
    if (ALPAQA_WITH_JSON)
        target_sources(ipopt-adapter PRIVATE
            "interop/ipopt/src/ipopt-json.cpp")
    endif()
    target_link_libraries(ipopt-adapter
        PUBLIC Ipopt::Ipopt alpaqa::alpaqa PRIVATE alpaqa::warnings)
    target_include_directories(ipopt-adapter PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/ipopt/include>
        $<INSTALL_INTERFACE:${ALPAQA_INSTALL_INCLUDEDIR}>)
    target_compile_definitions(ipopt-adapter PUBLIC ALPAQA_WITH_IPOPT)
    set_property(TARGET ipopt-adapter
        PROPERTY OUTPUT_NAME alpaqa-ipopt-adapter)
    alpaqa_configure_visibility(ipopt-adapter)
    set_target_properties(ipopt-adapter PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::ipopt-adapter ALIAS ipopt-adapter)
    list(APPEND ALPAQA_INSTALL_TARGETS ipopt-adapter)
    generate_export_header(ipopt-adapter
        EXPORT_FILE_NAME export/alpaqa/ipopt-adapter-export.h)
    target_include_directories(ipopt-adapter PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)
endif()

find_package(lbfgspp QUIET)
if (TARGET lbfgspp)
    message(STATUS "Compiling LBFGS++ solver adapter")
    add_library(lbfgspp-adapter "interop/lbfgspp/src/lbfgsb-adapter.cpp")
    target_link_libraries(lbfgspp-adapter
        PUBLIC lbfgspp alpaqa::alpaqa PRIVATE alpaqa::warnings)
    target_include_directories(lbfgspp-adapter PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/lbfgspp/include>
        $<INSTALL_INTERFACE:${ALPAQA_INSTALL_INCLUDEDIR}>)
    set_property(TARGET lbfgspp-adapter
        PROPERTY OUTPUT_NAME alpaqa-lbfgspp-adapter)
    alpaqa_configure_visibility(lbfgspp-adapter)
    set_target_properties(lbfgspp-adapter PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::lbfgspp-adapter ALIAS lbfgspp-adapter)
    list(APPEND ALPAQA_INSTALL_TARGETS lbfgspp-adapter)
    generate_export_header(lbfgspp-adapter
        EXPORT_FILE_NAME export/alpaqa/lbfgspp-adapter-export.h)
    target_include_directories(lbfgspp-adapter PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)
endif()

if (ALPAQA_WITH_LBFGSB)
    message(STATUS "Compiling L-BFGS-B solver adapter")
    add_library(lbfgsb-fortran STATIC
        "thirdparty/lbfgsb/Lbfgsb.3.0/lbfgsb.f"
        "thirdparty/lbfgsb/Lbfgsb.3.0/timer.f"
        "thirdparty/lbfgsb/Lbfgsb.3.0/blas.f"
        "thirdparty/lbfgsb/Lbfgsb.3.0/linpack.f"
    )
    set_target_properties(lbfgsb-fortran PROPERTIES LINKER_LANGUAGE Fortran)
    list(APPEND ALPAQA_INSTALL_TARGETS lbfgsb-fortran)
    add_library(lbfgsb-adapter
        "interop/lbfgsb/src/lbfgsb-adapter.cpp"
        "interop/lbfgsb/src/lbfgsb_c.F90"
    )
    if (ALPAQA_WITH_DRIVERS)
        target_sources(lbfgsb-adapter PRIVATE
            "interop/lbfgsb/src/lbfgsb-params.cpp")
    endif()
    if (ALPAQA_WITH_JSON)
        target_sources(lbfgsb-adapter PRIVATE
            "interop/lbfgsb/src/lbfgsb-json.cpp")
    endif()
    set_property(SOURCE "interop/lbfgsb/src/lbfgsb_c.F90"
                 PROPERTY Fortran_PREPROCESS On)
    target_link_libraries(lbfgsb-adapter
        PUBLIC lbfgsb-fortran alpaqa::alpaqa PRIVATE alpaqa::warnings)
    target_include_directories(lbfgsb-adapter PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/lbfgsb/include>
        $<INSTALL_INTERFACE:${ALPAQA_INSTALL_INCLUDEDIR}>)
    target_compile_definitions(lbfgsb-adapter PUBLIC ALPAQA_WITH_LBFGSB)
    set_property(TARGET lbfgsb-adapter
        PROPERTY OUTPUT_NAME alpaqa-lbfgsb-adapter)
    alpaqa_configure_visibility(lbfgsb-adapter)
    set_target_properties(lbfgsb-adapter PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::lbfgsb-adapter ALIAS lbfgsb-adapter)
    list(APPEND ALPAQA_INSTALL_TARGETS lbfgsb-adapter)
    generate_export_header(lbfgsb-adapter
        EXPORT_FILE_NAME export/alpaqa/lbfgsb-adapter-export.h)
    target_include_directories(lbfgsb-adapter PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)
endif()

if (ALPAQA_WITH_QPALM)
    message(STATUS "Compiling QPALM solver adapter")
    find_package(QPALM REQUIRED COMPONENTS cxx)
    add_library(qpalm-adapter "interop/qpalm/src/qpalm-adapter.cpp")
    if (ALPAQA_WITH_DRIVERS)
        target_sources(qpalm-adapter PRIVATE
            "interop/qpalm/src/qpalm-params.cpp")
    endif()
    if (ALPAQA_WITH_JSON)
        target_sources(qpalm-adapter PRIVATE
            "interop/qpalm/src/qpalm-json.cpp")
    endif()
    target_link_libraries(qpalm-adapter
        PUBLIC QPALM::qpalm_cxx alpaqa::alpaqa PRIVATE alpaqa::warnings)
    target_include_directories(qpalm-adapter PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interop/qpalm/include>
        $<INSTALL_INTERFACE:${ALPAQA_INSTALL_INCLUDEDIR}>)
    target_compile_definitions(qpalm-adapter PUBLIC ALPAQA_WITH_QPALM)
    set_property(TARGET qpalm-adapter
        PROPERTY OUTPUT_NAME alpaqa-qpalm-adapter)
    alpaqa_configure_visibility(qpalm-adapter)
    set_target_properties(qpalm-adapter PROPERTIES SOVERSION ${PROJECT_VERSION})
    add_library(alpaqa::qpalm-adapter ALIAS qpalm-adapter)
    list(APPEND ALPAQA_INSTALL_TARGETS qpalm-adapter)
    generate_export_header(qpalm-adapter
        EXPORT_FILE_NAME export/alpaqa/qpalm-adapter-export.h)
    target_include_directories(qpalm-adapter PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)
endif()

# Drivers
if (ALPAQA_WITH_DRIVERS)
    add_executable(driver
        "alpaqa/src/driver/alpaqa-driver.cpp"
        "alpaqa/src/driver/panoc-driver.cpp"
        "alpaqa/src/driver/fista-driver.cpp"
        "alpaqa/src/driver/pantr-driver.cpp"
        "alpaqa/src/driver/ipopt-driver.cpp"
        "alpaqa/src/driver/lbfgsb-driver.cpp"
        "alpaqa/src/driver/qpalm-driver.cpp"
        "alpaqa/src/driver/problem.cpp"
        "alpaqa/src/driver/param-complete.cpp"
        "alpaqa/src/driver/openmp.cpp"
    )
    set_target_properties(driver PROPERTIES
        OUTPUT_NAME "alpaqa-driver"
        RELEASE_POSTFIX ""
        DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}"
        RELWITHDEBINFO_POSTFIX ""
        MINSIZEREL_POSTFIX ""
    )
    target_link_libraries(driver
        PRIVATE alpaqa::alpaqa alpaqa::warnings)
    add_executable(alpaqa::driver ALIAS driver)
    list(APPEND ALPAQA_INSTALL_EXE driver)
    if (TARGET alpaqa::dl-loader)
        target_link_libraries(driver PRIVATE alpaqa::dl-loader)
    endif()
    if (TARGET alpaqa::casadi-loader)
        target_link_libraries(driver PRIVATE alpaqa::casadi-loader)
    endif()
    if (ALPAQA_WITH_EXTERNAL_CASADI)
        target_link_libraries(driver PRIVATE casadi::casadi)
    endif()
    if (TARGET alpaqa::cutest-interface)
        target_link_libraries(driver PRIVATE alpaqa::cutest-interface)
    endif()
    if (TARGET alpaqa::casadi-ocp-loader)
        target_link_libraries(driver PRIVATE alpaqa::casadi-ocp-loader)
    endif()
    if (TARGET alpaqa::ipopt-adapter)
        target_link_libraries(driver PRIVATE alpaqa::ipopt-adapter)
    endif()
    if (TARGET alpaqa::lbfgsb-adapter)
        target_link_libraries(driver PRIVATE alpaqa::lbfgsb-adapter)
    endif()
    if (TARGET alpaqa::qpalm-adapter)
        target_link_libraries(driver PRIVATE alpaqa::qpalm-adapter)
    endif()
endif()

# Gradient checker tool
if (ALPAQA_WITH_GRADIENT_CHECKER)
    add_executable(gradient-checker
        "alpaqa/src/driver/gradient-checker.cpp"
        "alpaqa/src/driver/problem.cpp"
    )
    set_target_properties(gradient-checker PROPERTIES
        OUTPUT_NAME "alpaqa-gradient-checker"
        RELEASE_POSTFIX ""
        DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}"
        RELWITHDEBINFO_POSTFIX ""
        MINSIZEREL_POSTFIX ""
    )
    target_link_libraries(gradient-checker
        PRIVATE alpaqa::alpaqa alpaqa::warnings)
    add_executable(alpaqa::gradient-checker ALIAS gradient-checker)
    list(APPEND ALPAQA_INSTALL_EXE gradient-checker)
    if (TARGET alpaqa::dl-loader)
        target_link_libraries(gradient-checker PRIVATE alpaqa::dl-loader)
    endif()
    if (TARGET alpaqa::casadi-loader)
        target_link_libraries(gradient-checker PRIVATE alpaqa::casadi-loader)
    endif()
    if (ALPAQA_WITH_EXTERNAL_CASADI)
        target_link_libraries(gradient-checker PRIVATE casadi::casadi)
    endif()
    if (TARGET alpaqa::cutest-interface)
        target_link_libraries(gradient-checker PRIVATE alpaqa::cutest-interface)
    endif()
    if (TARGET alpaqa::casadi-ocp-loader)
        target_link_libraries(gradient-checker PRIVATE alpaqa::casadi-ocp-loader)
    endif()
endif()

# Installation
include(cmake/Install.cmake)
