find_package(Eigen3 REQUIRED)

function(alpaqa_configure_visibility target)
    set_target_properties(${target} PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                                               VISIBILITY_INLINES_HIDDEN true)
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
        target_link_options(${target} PRIVATE "LINKER:--exclude-libs,ALL")
    endif()
endfunction()

add_library(alpaqa
    "src/util/type-erasure.cpp"
    "src/util/demangled-typename.cpp"
    "src/util/print.cpp"
    "src/accelerators/lbfgs.cpp"
    "src/problem/problem-counters.cpp"
    "src/problem/ocproblem-counters.cpp"
    "src/outer/alm.cpp"
    "src/inner/panoc.cpp"
    "src/inner/panoc-ocp.cpp"
    "src/panoc-alm.cpp"
    "src/structured-panoc-alm.cpp"
    "src/inner/internal/solverstatus.cpp"
    "src/inner/directions/panoc/structured-lbfgs.cpp"
    "src/inner/internal/panoc-helpers.cpp"
    "include/alpaqa/util/type-erasure.hpp"
)
target_compile_features(alpaqa PUBLIC cxx_std_20)
target_include_directories(alpaqa PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_definitions(alpaqa PUBLIC
    $<$<CONFIG:Debug>:EIGEN_INITIALIZE_MATRICES_BY_NAN>
    $<$<CONFIG:Debug>:EIGEN_RUNTIME_NO_MALLOC>)
target_link_libraries(alpaqa PUBLIC Eigen3::Eigen)
target_link_libraries(alpaqa PRIVATE warnings)
alpaqa_configure_visibility(alpaqa)
target_compile_definitions(alpaqa PUBLIC
    $<$<BOOL:${ALPAQA_WITH_QUAD_PRECISION}>:ALPAQA_WITH_QUAD_PRECISION>)
target_link_libraries(alpaqa PUBLIC
    $<$<BOOL:${ALPAQA_WITH_QUAD_PRECISION}>:quadmath>)
add_library(alpaqa::alpaqa ALIAS alpaqa)
list(APPEND ALPAQA_INSTALL_TARGETS alpaqa)

# DLL import/export
include(GenerateExportHeader)
generate_export_header(alpaqa
    EXPORT_FILE_NAME export/alpaqa/export.h)
target_include_directories(alpaqa PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)

# CasADi
if (ALPAQA_WITH_CASADI)
    find_package(casadi REQUIRED)
    find_package(Threads REQUIRED)
    # Normal NLPs
    add_library(casadi-loader "interop/casadi/CasADiProblem.cpp")
    set_property(TARGET casadi-loader
        PROPERTY OUTPUT_NAME alpaqa-casadi-loader)
    target_link_libraries(casadi-loader
        PUBLIC  alpaqa::alpaqa
        PRIVATE casadi Threads::Threads alpaqa::warnings)
    alpaqa_configure_visibility(casadi-loader)
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

set(ALPAQA_WITH_CASADI_OCP_SUPPORTED true)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND
    CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16.0.0)
    set(ALPAQA_WITH_CASADI_OCP_SUPPORTED false)
endif()

include(CMakeDependentOption)
cmake_dependent_option(ALPAQA_WITH_CASADI_OCP
    "Build the CasADi loader for optimal control problems"
    ${ALPAQA_WITH_CASADI_OCP_SUPPORTED} "ALPAQA_WITH_CASADI" OFF)

# Optimal control problems
if (ALPAQA_WITH_CASADI_OCP)
    add_library(casadi-ocp-loader "interop/casadi/CasADiControlProblem.cpp")
    set_property(TARGET casadi-ocp-loader
        PROPERTY OUTPUT_NAME alpaqa-casadi-ocp-loader)
    target_link_libraries(casadi-ocp-loader
        PUBLIC  alpaqa::alpaqa
        PRIVATE casadi Threads::Threads alpaqa::warnings)
    alpaqa_configure_visibility(casadi-ocp-loader)
    add_library(alpaqa::casadi-ocp-loader ALIAS casadi-ocp-loader)
    list(APPEND ALPAQA_INSTALL_TARGETS casadi-ocp-loader)
    generate_export_header(casadi-ocp-loader
        EXPORT_FILE_NAME export/alpaqa/casadi-ocp-loader-export.h)
    target_include_directories(casadi-ocp-loader PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)
endif()

if (CMAKE_DL_LIBS)
    add_library(dl-headers INTERFACE)
    # TODO: move to independent folder
    target_include_directories(dl-headers INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
    add_library(alpaqa::dl-headers ALIAS dl-headers)
    list(APPEND ALPAQA_INSTALL_TARGETS dl-headers)
    add_library(dl-loader "interop/dl/dl-problem.cpp")
    set_property(TARGET dl-loader
        PROPERTY OUTPUT_NAME alpaqa-dl-loader)
    target_link_libraries(dl-loader
        PUBLIC alpaqa::alpaqa alpaqa::dl-headers
        PRIVATE ${CMAKE_DL_LIBS} alpaqa::warnings)
    add_library(alpaqa::dl-loader ALIAS dl-loader)
    list(APPEND ALPAQA_INSTALL_TARGETS dl-loader)
    message(STATUS "Compiling WITH dl support")
else()
    message(STATUS "Compiling WITHOUT dl support")
endif()

# Installation
include(cmake/Install.cmake)
